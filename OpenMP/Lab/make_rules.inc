#
# Makefile rules
#

CC = gcc
FC = gfortran
# for intel using icc and ifort
#CC = icc
#FC = ifort

# Setup debug and optimized versions of the code

OPTDIR = Optimized
DEBUGDIR = Debug

OPENMPFLAGS = -fopenmp
ifeq "$(findstring release, $(MAKECMDGOALS))" ""
	BUILD=$(DEBUGDIR)
	CCFLAGS += -g3 -O0 -Wall -Wno-unknown-pragmas -Wno-format $(OPENMPFLAGS)
	FCFLAGS += -g3 -O0 -Wall -Wno-unknown-pragmas -Wno-format $(OPENMPFLAGS)
else
	BUILD=$(OPTDIR)
	CCFLAGS += -O3 $(OPENMPFLAGS)
	FCFLAGS += -O3 $(OPENMPFLAGS)
endif


# Note:
# The opt-report levels for Intel compilers are:
# 0: No vectorization report
# 1: Reports when vectorization has occurred.
# 2: Adds diagnostics why vectorization did not occur.
# 3: Adds vectorization loop summary diagnostics.
# 4: Additional detail, e.g. on data alignment
# 5: Adds detailed data dependency information

.SUFFIXES:
.SUFFIXES: .c .f90 .o .mod
.PHONEY:   clean


$(BUILD)/%.o: %.f90 Makefile ../make_rules.inc
	$(FC) $(FCFLAGS) -c -o $@ $<

$(BUILD)/%.o: %.c Makefile ../make_rules.inc
	$(CC) $(CCFLAGS) -c -o $@ $<


all: $(EXES)

clean:
	$(RM) *.o *.mod $(EXES)
