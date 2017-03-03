#
# Makefile rules
#

CC      = mpicc
FC      = mpifort
DEBUG   = -g -O0
CFLAGS  = -std=c11 $(DEBUG)
FCFLAGS = $(DEBUG)


# Note:

.SUFFIXES:
.SUFFIXES: .c .f90 .o .mod
.PHONEY:   clean

%:: %.f90
	$(FC) -o $@ $(FCFLAGS) $<

%:: %.F90
	$(FC) -fpp -o $@ $(FCFLAGS) $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

all: $(EXES)

clean:
	$(RM) *.o *.mod $(EXES)
