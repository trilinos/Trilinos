# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBPETRA_TSF) by typing 'make'.
#


TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

LIBPETRA_TSF = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libpetra_tsf.a

include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

# Petra communication defines
PETRA_COMM_SERIAL          = SERIAL
PETRA_COMM_MPI             = PETRA_MPI
PETRA_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))

# Petra communication defines
TSF_COMM_SERIAL          = TSF_SERIAL
TSF_COMM_MPI             = TSF_MPI
TSF_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM) -D$(TSF_COMM)

PETRA_INCLUDES = -I../petra
TSF_INCLUDES = -I../framework
INCLUDES = -I.$(ARCH_INCLUDES) $(PETRA_INCLUDES) $(TSF_INCLUDES) $(BLAS_INCLUDES) 

LIB_PATHS=

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
#=======================================================================
# PETRA source files
#=======================================================================


PETRA_TSF_CC = Petra_TSF_RDP_MultiVector.cc
PETRA_TSF_C  =
PETRA_TSF_F =

#=======================================================================
# PETRA include files
#=======================================================================

PETRA_TSF_INC =

PETRA_TSF_OBJ          = $(PETRA_TSF_CC:.cc=.o) $(PETRA_TSF_C:.c=.o) $(PETRA_TSF_F:.f=.o)

$(LIBPETRA_TSF): $(PETRA_TSF_OBJ)
	mkdir -p $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)
	$(AR) rcv $(LIBPETRA_TSF) $(PETRA_TSF_OBJ)
	$(RANLIB) $(LIBPETRA_TSF)

#
#include ../../etc/depends.petra_tsf

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  $(LIBPETRA) *~
