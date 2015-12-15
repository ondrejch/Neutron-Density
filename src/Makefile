# Simple makefile
# Ondrej Chvala, ochvala@utk.edu, 20151212
# GNU/GPL
#
EXECUTABLE=../bin/main.exe
SOURCES=inputinterpreter.f90  neutrondens.f90 main.f90

OBJECTS=$(SOURCES:.f90=.o)
LDFLAGS=-llapack

FC=gfortran
CFLAGS=-g -finit-local-zero -fcheck=all -ftrapv -ffpe-trap=overflow,invalid -Wall -Wextra
#CFLAGS=-O3 -Wall -Wextra 

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(FC) $(OBJECTS) $(LDFLAGS) -o $@

%.o: %.f90
	$(FC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.mod $(EXECUTABLE)

