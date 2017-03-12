PROG = ex3
include ../../../build_flags.inc



$(PROG):$(PROG).f90
	$(FC) $(FCFLAGS) -c timing.f90
	$(FC) $(FCFLAGS) -o $(PROG).out $< timing.o
	$(FC) $(DBGFLAGS) -o $(PROG).dbg $< timing.o
all: $(PROG)
clean:
	@rm -rf *.out *.dbg
	@rm -rf *.o *.mod
