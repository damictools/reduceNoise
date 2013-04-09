#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>

#include "globalConstants.h"


#include "TFile.h"
#include "TNtuple.h"
#include "TPrincipal.h"

using namespace std;

const int kVerbose = 1;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}

void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program computes overscan mean and subtracts it line by line.\n";
    cout << "It handles all the available HDUs. The HDUs in the output fit file\n";
    cout << "will be double (64bits) for 64bits images and float (32bis) in all\n";
    cout << "the other cases.\n";
    cout << "The card \"TRIMSEC\" must be present in the header to use this program.\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file> -o <output filename> \n\n";
  cout << "\nOptions:\n";
  cout << "  -v for verbosity\n";
  cout << "  -s <HDU number> for processing a single HDU \n\n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

string bitpix2TypeName(int bitpix){
  
  string typeName;
  switch(bitpix) {
      case BYTE_IMG:
          typeName = "BYTE(8 bits)";
          break;
      case SHORT_IMG:
          typeName = "SHORT(16 bits)";
          break;
      case LONG_IMG:
          typeName = "INT(32 bits)";
          break;
      case FLOAT_IMG:
          typeName = "FLOAT(32 bits)";
          break;
      case DOUBLE_IMG:
          typeName = "DOUBLE(64 bits)";
          break;
      default:
          typeName = "UNKNOWN";
  }
  return typeName;
}

void computeMVASigma(const int iExt, const float cutVal, const vector< double* > &vPix, const vector< int > &vFullNCol, const vector< int > &vNLines, vector<double> &vCoef){
  
  const unsigned int kN = vPix.size();
  Float_t vR[kN];
  Float_t vL[kN];
  
  TPrincipal* principal = new TPrincipal(kN+1,"");
  
  const long fullNCol = vFullNCol[0];
  const long nCol     = vFullNCol[0]/2;
  const long nLines   = vNLines[0];
  const Long64_t nPix = nCol*nLines;

  Long64_t npri=0;
  if(kVerbose){
    cout << "Reading data..\n";
  }
  
  for(int l=0;l<nLines;++l){
    for(int c=0;c<nCol;++c){
      
      Double_t data[kN+1];
      
      data[kN] = vPix[iExt][fullNCol*l + c + nCol]; //vR[iExt-1];
      for(int j=0;j<kN;++j) data[j] = vPix[j][ fullNCol*l - c + nCol -1];
      
      bool skip=false;
      for(int j=0;j<kN+1;++j){
	if(fabs(data[j])>cutVal){
	  skip=true;
	  break;
	}
      }
      
      if(skip) continue;
      
      principal->AddRow(data);
      ++npri;
    }
  }
  
  if(kVerbose){
    cout << "Done\nComputing Covariance Matrix..";
  }
  const TMatrixD *m = (principal->GetCovarianceMatrix());
  
//   m->Print();
  TMatrixD MCovY(kN,kN);
  
  TMatrixD rCovXiY(kN,1);
  
  for(int i=0;i<kN;++i){
    rCovXiY[i][0]=(*m)[kN][i];
    for(int j=0;j<=i;++j){
       MCovY[i][j]= (*m)[i][j];
       
       if(j!=i) MCovY[j][i]+= (*m)[i][j];
    }
  }
  
  if(kVerbose){
    MCovY.Print();
    cout << "Done\nInverting Matrix..\n";
  }
  TMatrixD invMCovY(MCovY);
  invMCovY.Invert();
  
  TMatrixD ai(invMCovY*rCovXiY);

  ostringstream transForm;
  transForm << "R" <<iExt << " - (";
  for(int i=0;i<kN;++i){
    transForm << "+"<< ai[i][0] << "*L" << i;
    vCoef.push_back(ai[i][0]);
  }
  transForm << ")"; 
  
  cout << transForm.str() << endl;
  
}




int computeImage(const char *inFile, const char *outFile, const vector<int> &outHdu){
  int status = 0;
  double nulval = 0.;
  int anynul = 0;
  int  nhdu = 0;
  
  fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
  fits_open_file(&infptr, inFile, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  fits_get_num_hdus(infptr, &nhdu, &status);
  if (status != 0) return(status);
  
  vector<int> inHdu;
//   for(unsigned int i=0;i<inHdu.size();++i){
//     if(inHdu[i] > nhdu){
//       fits_close_file(infptr,  &status);
//       cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
//       return -1000;
//     }
//   }
  
  if(inHdu.size() == 0){
    for(int i=0;i<nhdu;++i){
      inHdu.push_back(i+1);
    }
  }
  const unsigned int nUseHdu=inHdu.size();
  
  
  fitsfile *outfptr;   /* FITS file pointers defined in fitsio.h */
  fits_create_file(&outfptr, outFile, &status);
  
  
  vector< double* > vPix;
  vector< int > vFullNCol;
  vector< int > vNLines;
  
  for(unsigned int eI=0; eI<nUseHdu; ++eI)  /* Main loop through each extension */
  {
    const unsigned int n = inHdu[eI];
    
    const int nHDUsToProcess = nUseHdu;
    
    /* get input image dimensions and total number of pixels in image */
    int hdutype, bitpix, bytepix, naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
    long totpix = naxes[0] * naxes[1];
    double bzero;
    ffgky(infptr, TDOUBLE, "BZERO", &bzero, NULL, &status);
    if (status){
      status = 0;
      bzero = 0.0;
    }
    
    /* Don't try to process data if the hdu is empty */    
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      continue;
    }
    
    
    double* inArray = new double[totpix];
    vPix.push_back(inArray);
    
    for(int i=0;i<totpix;++i) inArray[i] = 0;
    
    if(gVerbosity){
      showProgress((n-1)*3+0,3*nHDUsToProcess);
    }
      
    int xMin=1;
    int xMax=naxes[0];
    int yMin=1;
    int yMax=naxes[1];
    
    /* Read the images as doubles, regardless of actual datatype. */
    long fpixel[2]={xMin,yMin};
    long lpixel[2]={xMax,yMax};
    long inc[2]={1,1};
    fits_read_subset(infptr, TDOUBLE, fpixel, lpixel, inc, &nulval, inArray, &anynul, &status);
    if (status != 0) return(status);
    if(gVerbosity){
      showProgress((n-1)*3+1,3*nHDUsToProcess);
    }
    
    vFullNCol.push_back(naxes[0]);
    vNLines.push_back(naxes[1]);
    
    if(gVerbosity){
      showProgress((n-1)*3+2,3*nHDUsToProcess);
    }

  }
  
  
  

  const int nOutHdu = outHdu.size();
  vector<float*> vPixOut;
  for (unsigned int oi=0; oi<nOutHdu; ++oi)  /* Main loop through each extension */
  {
    const int i = outHdu[oi]-1;
    float cutVal = 200;
    vector<double> vCoef;
    computeMVASigma(i, cutVal, vPix, vFullNCol, vNLines, vCoef);
    
    /* Explicitly create new image */
    long naxesOut[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    naxesOut[0] = vFullNCol[i]/2;
    naxesOut[1] = vNLines[i];
    int bitpix = FLOAT_IMG;
    int naxis = 2;
    
    const long fullNCol = vFullNCol[0];
    const long nCol     = naxesOut[0];
    const long nLines   = naxesOut[1];
    const long totpixOut = naxesOut[0]*naxesOut[1];
    float* outArray = new float[totpixOut];
    vPixOut.push_back(outArray);
    for(int l=0;l<nLines;++l){
      for(int c=0;c<nCol;++c){
	outArray[nCol*l + c] = vPix[i][ fullNCol*l + c + nCol];// - vPix[i][ fullNCol*l - c + nCol -1]
	for(unsigned int ex=0; ex<vCoef.size(); ++ex){
	  outArray[nCol*l + c] -= vCoef[ex]*vPix[ex][ fullNCol*l - c + nCol -1];
	}
      }
    }
    
    
    
    /* Explicitly create new image */
    fits_create_img(outfptr, bitpix, naxis, naxesOut, &status);
    if (status) {
      fits_report_error(stderr, status);
      return(status);
    }
    /* copy the relevant keywords (not the structural keywords) */
    int hdutype;
    fits_movabs_hdu(infptr, outHdu[oi], &hdutype, &status);
    int nkeys = 0;
    fits_get_hdrspace(infptr, &nkeys, NULL, &status); 
    for (int r = 1; r <= nkeys; ++r) {
      char card[FLEN_CARD];
      fits_read_record(infptr, r, card, &status);
      if (fits_get_keyclass(card) > TYP_CMPRS_KEY) fits_write_record(outfptr, card, &status);
    }
    
    
    long first = 1;
    fits_write_img(outfptr, TFLOAT, first, totpixOut, outArray, &status);
  }
   
  
  /* Close the fits files */
  fits_close_file(infptr,  &status);
  
  fits_close_file(outfptr,  &status);
  
  /* clean up */
  for(unsigned int i=0;i<vPix.size();++i){
    delete[] vPix[i];
  }
  for(unsigned int i=0;i<vPixOut.size();++i){
    delete[] vPixOut[i];
  }
  
  if(gVerbosity){
    showProgress(1,1);
  }
  
  return status;
}


void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], vector<int> &singleHdu, string &inFile, string &outFile){
  
  if(argc == 1) return 1;
  
  bool outFileFlag = false;
  int opt=0;
  singleHdu.clear();
  while ( (opt = getopt(argc, argv, "i:o:s:vVhH?")) != -1) {
    switch (opt) {
    case 'o':
      if(!outFileFlag){
        outFile = optarg;
        outFileFlag = true;
      }
      else{
        cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
        return 2;
      }
      break;
    case 's':
      singleHdu.push_back(atoi(optarg));
      break;
    case 'V':
    case 'v':
      gVerbosity = 1;
      break;
    case 'h':
    case 'H':
    default: /* '?' */
      return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }

  inFile="";
  
  if(argc-optind==0){
    cerr << red << "Error: no input file(s) provided!\n\n" << normal;
    return 1;
  }
  else if(argc-optind>1){
    cerr << red << "Error: more than one input file provided!\n\n" << normal;
    return 1;
  }
  
  inFile=argv[optind];
  if(!fileExist(inFile.c_str())){
    cout << red << "\nError reading input file: " << inFile <<"\nThe file doesn't exist!\n\n" << normal;
    return 1;
  }
  
  std::sort(singleHdu.begin(), singleHdu.end());
  
  return 0;
}

int main(int argc, char *argv[])
{
  
  checkArch(); //Check the size of the double and float variables.
  
  time_t start,end;
  double dif;
  time (&start);
  
  string outFile;
  string inFile;
  vector<int> singleHdu;
  
  int returnCode = processCommandLineArgs( argc, argv, singleHdu, inFile, outFile);
  if(returnCode!=0){
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }
  
  /* Overwrite the output file if it already exist */
  if(fileExist(outFile.c_str())){
    cout << yellow << "\nThe output file exist. " << normal;
    deleteFile(outFile.c_str());
  }
  
  if(gVerbosity){
    cout << bold << "\nWill read the following file:\n" << normal;
    cout << "\t" << inFile << endl;
    cout << bold << "\nThe output will be saved in the file:\n\t" << normal << outFile << endl;
  }
  
  int status = computeImage( inFile.c_str(),  outFile.c_str(), singleHdu);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) cout << green << "All done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  return status;
}
