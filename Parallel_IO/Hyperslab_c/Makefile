#
# MPI IO test programs
#

CFLAGS   = -O3
F90FLAGS = -O3

CC      = h5pcc
HFC      = h5pfc

EXES = hyperslab_by_row read_hyperslab

.SUFFIXES:
.SUFFIXES:    .c .f90 .o
.PHONEY:      clean

all: $(EXES)

hyperslab_by_row: hyperslab_by_row.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

read_hyperslab: read_hyperslab.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


clean:
	-$(RM) *.h5 *.o *.mod $(EXES)

