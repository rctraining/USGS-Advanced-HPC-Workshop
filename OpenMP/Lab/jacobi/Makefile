#
# Makefile for OpenMP Jacobian Solver
#

FC      = ifort
CC      = icc
DEBUG   = -g -O3
REPORT  = -qopt-report \
          -qopt-report-phase=all
          -qopt-report-phase=vec,loop
FCFLAGS = -fopenmp $(REPORT) $(DEBUG) \
          -D__BIGGEST_ALIGNMENT__=64  \

EXE  = jacobi
SRCS = kinds.f90     \
       err.f90       \
       domain.f90    \
       args.f90      \
       timing.c      \
       walltime.f90  \
       jacobi.f90    \
       main.f90
#       hdf_rw.f90
OBJS = kinds.o       \
       err.o         \
       domain.o      \
       args.o        \
       timing.o      \
       walltime.o    \
       jacobi.o      \
       main.o
#       hdf_rw.o
#OBJS := $(patsubst %.f90,%.o,$(SRCS)) \
#        $(patsubst %.c,%.o,$(SRCS))

.SUFFIXES:
.SUFFIXES: .c .f90 .d .o

$(EXE): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $(OBJS) $(LDFLAGS)

clean:
	-$(RM) $(EXE) *.o *.mod

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
