PROG1 = ex4a
PROG2 = ex4b
include ../../../build_flags.inc

all: $(PROG1) $(PROG2)


$(PROG1):$(PROG1).f90
	$(FC) $(FCFLAGS) -c timing.f90
	$(FC) $(FCFLAGS) -o $(PROG1).out $< timing.o
	$(FC) $(DBGFLAGS) -o $(PROG1).dbg $< timing.o

$(PROG2):$(PROG2).f90
	$(FC) $(FCFLAGS) -c timing.f90
	$(FC) $(FCFLAGS) -o $(PROG2).out $< timing.o
	$(FC) $(DBGFLAGS) -o $(PROG2).dbg $< timing.o

clean:
	@rm -rf *.out *.dbg
	@rm -rf *.o *.mod
