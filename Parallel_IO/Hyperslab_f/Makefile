#
# MPI IO test programs
#

CFLAGS   = -O3
F90FLAGS = -O3

HCC      = h5pcc
HFC      = h5pfc

EXES = hyperslab_by_col read_hyperslab

.SUFFIXES:
.SUFFIXES:    .c .f90 .o
.PHONEY:      clean

all: $(EXES)

hyperslab_by_col: hyperslab_by_col.f90
	$(HFC) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

read_hyperslab: read_hyperslab.f90
	$(HFC) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

.f90.o:
	-$(RM) $*.o $*.mod
	$(HFC) -c $(F90FLAGS) -o $@ $<

clean:
	-$(RM) *.h5 *.o *.mod $(EXES)

