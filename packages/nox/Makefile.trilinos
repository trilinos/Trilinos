# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBNLSPACK) by typing 'make'.
#
include $(TRILINOS_HOME)/build/TARGET

TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

LIBNLSPACK = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libnlspack.a

include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

# Petra communication defines
PETRA_COMM_SERIAL          = SERIAL
PETRA_COMM_MPI             = PETRA_MPI
PETRA_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))

# NLSPACK communication defines
NLSPACK_COMM_SERIAL          = SERIAL
NLSPACK_COMM_MPI             = AZTEC_MPI
NLSPACK_COMM                 = $(NLSPACK_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(NLSPACK_ARCH_DEFINES) -D$(NLSPACK_COMM) \
         -D$(PETRA_COMM) $(PETRA_ARCH_DEFINES) -DIFPACK -DSPARSEBLAS \
         -DAZ_COL_REORDER

INCLUDES = -I. $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src \
           -I$(TRILINOS_HOME)/packages/aztecoo/src

LIB_PATHS=

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS = $(CXXFLAGS)
#=======================================================================
# NLSPACK source files
#=======================================================================

NLSPACK_CC = NLS_Parameter.C NLS_ParameterList.C \
             NLS_PetraVector.C NLS_PetraGroup.C \
             NLS_MethodManager.C NLS_Newton.C

NLSPACK_C  = 

#=======================================================================
# NLSPACK include files
#=======================================================================

NLSPACK_INC =

NLSPACK_OBJ          = $(NLSPACK_CC:.C=.o)  $(NLSPACK_C:.c=.o)

$(LIBNLSPACK): $(NLSPACK_OBJ)
	mkdir -p $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)
	$(AR) rcv $(LIBNLSPACK) $(NLSPACK_OBJ)
	$(RANLIB) $(LIBNLSPACK)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.ifpack

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  $(LIBNLSPACK) *~
