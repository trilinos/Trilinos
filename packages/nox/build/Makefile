# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBNLSPACK) by typing 'make'.
#

include Makefile.defs

all: $(LIBNLSPACK)
	$(DO_NOTHING)

#=======================================================================
# NLSPACK source files
#=======================================================================

include Makefile.NLSPACK.autorules

$(LIBNLSPACK): $(NLSPACK_OBJ)
	mkdir -p $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)
	$(AR) rcv $(LIBNLSPACK) $(NLSPACK_OBJ)
	$(RANLIB) $(LIBNLSPACK)

clean:
	@rm -f *.o  $(LIBNLSPACK) *~
