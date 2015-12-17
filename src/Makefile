# Simple makefile
# Ondrej Chvala, ochvala@utk.edu, 20151212
# GNU/GPL
#
EXEDIR     = ../bin
EXECUTABLE = $(EXEDIR)//main.exe
SOURCES    = feedback.f90 inputinterpreter.f90  neutrondens.f90 main.f90

OBJECTS    = $(SOURCES:.f90=.o)
LDFLAGS    = -llapack

FC         = gfortran
CFLAGS     = -g -finit-local-zero -fcheck=all -ftrapv -ffpe-trap=overflow,invalid -Wall -Wextra
#CFLAGS     = -O3 -Wall -Wextra 

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	@if [ ! -d $(EXEDIR) ]; then mkdir -p $(EXEDIR); fi
	$(FC) $(OBJECTS) $(LDFLAGS) -o $@

%.o: %.f90
	$(FC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.mod $(EXECUTABLE)

