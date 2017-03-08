PROG = ex5
GPU_OUT = gpu
SERIAL_OUT = serial
FC = pgfortran
FCFLAGS_GPU = -acc -fast -ta=tesla -Minfo=accel -mp
FCFLAGS_SERIAL = -fast -mp

$(PROG): $(PROG).f90
	$(FC) $(FCFLAGS_GPU) -o $(PROG).$(GPU_OUT) $<
	$(FC) $(FCFLAGS_SERIAL) -o $(PROG).$(SERIAL_OUT) $<
all: $(PROG)
clean:
	@rm -rf *.gpu *.serial
