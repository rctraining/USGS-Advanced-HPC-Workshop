SERIAL = ex1
FAUX = fparallel
GPU0 = gpu0
GPU1 = gpu1
GPU2 = gpu2
GPU3 = gpu3
OUT = out

FC = pgfortran
FCFLAGS_SERIAL = -fast
FCFLAGS_GPU = -fast -ta=tesla

ex1: ex1.f90
	$(FC) $(FCFLAGS_SERIAL) -o $(SERIAL).$(OUT) $<
fparallel: fparallel.f90
	$(FC) $(FCFLAGS_SERIAL) -o $(FAUX).$(OUT) $<

gpu0:   gpu0.f90
	$(FC) $(FCFLAGS_GPU) -o $(GPU0).$(OUT) $<
gpu1:   gpu1.f90
	$(FC) $(FCFLAGS_GPU) -o $(GPU1).$(OUT) $<
gpu2:   gpu2.f90
	$(FC) $(FCFLAGS_GPU) -o $(GPU2).$(OUT) $<
gpu3:   gpu3.f90
	$(FC) $(FCFLAGS_GPU) -o $(GPU3).$(OUT) $<
all: clean serial fparallel gpu0 gpu1 gpu2 gpu3
clean:
	@rm -rf *.out
