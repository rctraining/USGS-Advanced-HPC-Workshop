PROG = ex4
GPU_OUT = gpu
SERIAL_OUT = serial
FC = pgfortran
FCFLAGS_GPU = -acc -fast -ta=tesla -Minfo=accel
FCFLAGS_SERIAL = -fast

$(PROG): $(PROG).f90
	$(FC) $(FCFLAGS_GPU) -o $(PROG).$(GPU_OUT) $<
	$(FC) $(FCFLAGS_SERIAL) -o $(PROG).$(SERIAL_OUT) $<
all: $(PROG)
clean:
	@rm -rf *.gpu *.serial
