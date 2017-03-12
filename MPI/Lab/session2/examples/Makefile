PROG = ex1
include ../../build_flags.inc

all: messages_blocking  messages_nonblocking

messages_blocking: messages_blocking.f90
	$(FC) $(FCFLAGS) -o messages_blocking.out $<
	$(FC) $(DBGFLAGS) -o messages_blocking.dbg $<

messages_nonblocking: messages_nonblocking.f90
	$(FC) $(FCFLAGS) -o messages_nonblocking.out $<
	$(FC) $(DBGFLAGS) -o messages_nonblocking.dbg $<

clean:
	@rm -rf *.out *.dbg
