# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or 
#    pclinux.
#       Other machines require an appropriate 
#       makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same
#       architecture and communication mode.
#
#    4) Make the archive $(LIBSPBLAS) by typing 'make'.
#

include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH)

INCLUDES =
LIB_PATHS=

CFLAGS=$(ARCH_CFLAGS) $(AZTEC_ARCH_DEFINES) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)

#=======================================================================
# SPBLAS source files
#=======================================================================

SPBLAS_C = cblas_duscr_vbr.c cblas_dusmm.c

#=======================================================================
# SPBLAS include files
#=======================================================================

SPBLAS_INC = spblas.h

SPBLAS_OBJ          = $(SPBLAS_C:.c=.$(OBJ_EXT))

#$(LIBSPBLAS): $(SPBLAS_OBJ)
LIBSPBLAS: $(SPBLAS_OBJ)
	mkdir -p $(TRILINOS_HOME)/external/spblas/lib/$(TRILINOS_TARGET)
	$(AR)  $(ARFLAGS) $(LIB_WIN_FLAGS)$(LIBSPBLAS) $(SPBLAS_OBJ)
# following ranlib commented out by rbl
#	$(RANLIB) $(LIBSPBLAS)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.spblas

clean:
	echo "junk" > dummy.$(OBJ_EXT)
	rm -f *.$(OBJ_EXT)

clobber:
	echo "junk" > dummy.$(OBJ_EXT)
	rm -f *.$(OBJ_EXT)  $(LIBSPBLAS)
