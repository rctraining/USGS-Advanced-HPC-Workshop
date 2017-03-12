PROG = ex1
FC = mpif90
FCFLAGS = -O3 -mavx -ffixed-line-length-132
DBGFLAGS = -O0 -fbounds-check -fbacktrace -ffixed-line-length-132

all: messages_blocking  messages_nonblocking

messages_blocking: messages_blocking.f90
	$(FC) $(FCFLAGS) -o messages_blocking.out $<
	$(FC) $(DBGFLAGS) -o messages_blocking.dbg $<

messages_nonblocking: messages_nonblocking.f90
	$(FC) $(FCFLAGS) -o messages_nonblocking.out $<
	$(FC) $(DBGFLAGS) -o messages_nonblocking.dbg $<

clean:
	@rm -rf *.out *.dbg
