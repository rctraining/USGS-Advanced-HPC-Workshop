PROG = ex4
FC = mpif90
FCFLAGS = -O3 -mavx -ffixed-line-length-132
DBGFLAGS = -O0 -fbounds-check -fbacktrace -ffixed-line-length-132

$(PROG): $(PROG).f90
	$(FC) $(FCFLAGS) -o $(PROG).out $<
	$(FC) $(DBGFLAGS) -o $(PROG).dbg $<
all: $(PROG)
clean:
	@rm -rf *.out *.dbg
