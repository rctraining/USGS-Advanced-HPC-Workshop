PROG = ex1_sol
OUT = out
FC = pgfortran
FCFLAGS_GPU = -acc -fast -ta=tesla -Minfo=accel

$(PROG): $(PROG).f90
	$(FC) $(FCFLAGS_GPU) -o $(PROG).$(OUT) $<
all: $(PROG)
clean:
	@rm -rf *.out
