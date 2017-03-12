PROG = ex2
include ../../../build_flags.inc



$(PROG): $(PROG).f90
	$(FC) $(FCFLAGS) -o $(PROG).out $<
	$(FC) $(DBGFLAGS) -o $(PROG).dbg $<
all: $(PROG)
clean:
	@rm -rf *.out *.dbg
