#
# MPI IO test programs
#

CFLAGS   = -O3
F90FLAGS = -O3 #-align all \
#	   -fp-model strict -fpe0 -openmp

HCC      = h5pcc
HFC      = h5pfc

EXES = hdf_pwrite

.SUFFIXES:
.SUFFIXES:    .c .f90 .o
.PHONEY:      clean

all: $(EXES)

hdf_pwrite: kinds.o hdf_pwrite.o
	$(HFC) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

.f90.o:
	-$(RM) $*.o $*.mod
	$(HFC) -c $(F90FLAGS) -o $@ $<

clean:
	-$(RM) *.h5 *.o *.mod $(EXES)

