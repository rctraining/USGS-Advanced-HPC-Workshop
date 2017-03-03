# Toplevel makefile to recursively make all examples.

SUBDIRS = hello_c \
          hello_f \
	  ring_allreduce \
	  ring_blocking_correct \
	  ring_blocking_deadlock \
	  ring_nonblocking

.PHONY: all clean

all clean:
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir $@;    \
	done

