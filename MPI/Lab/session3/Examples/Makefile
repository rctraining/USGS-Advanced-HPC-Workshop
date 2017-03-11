PROG = ex1
FC = mpif90
FCFLAGS = -O3 -mavx -ffixed-line-length-132
DBGFLAGS = -O0 -fbounds-check -fbacktrace -ffixed-line-length-132


row_read: row_read.f90
	$(FC) $(FCFLAGS) -o row_read.out $<
	$(FC) $(DBGFLAGS) -o row_read.dbg $<
all: $(PROG)
clean:
	@rm -rf *.out *.dbg
