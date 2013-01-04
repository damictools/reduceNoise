CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS += 
OBJECTS = reduceNoise.o 
HEADERS = globalConstants.h

ALL : reduceNoise.exe
	@echo "Listo!"

reduceNoise.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o reduceNoise.exe $(LIBS) $(GLIBS) $(CFLAGS)

reduceNoise.o : reduceNoise.cc $(HEADERS)
	$(CPP) -c reduceNoise.cc -o reduceNoise.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
