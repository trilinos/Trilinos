# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBAZTECOO) by typing 'make'.
#

include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

LIBEPETRA= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libepetra.a
LIBAZTECOO= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libaztecoo.a
LIBIFPACK= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libifpack.a
LIBTRILINOS_UTIL= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libtriutils.a

include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

# Epetra communication defines
EPETRA_COMM_SERIAL          = SERIAL
EPETRA_COMM_MPI             = EPETRA_MPI
EPETRA_COMM                 = $(EPETRA_COMM_$(TRILINOS_COMM))

FORMAT=VBR

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM) \
         -DIFPACK -D$(FORMAT)

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src \
           -I$(TRILINOS_HOME)/packages/aztecoo/src \
           -I$(TRILINOS_HOME)/packages/triutils/src

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)



LIB_PATHS= $(LIBAZTECOO) $(LIBEPETRA) $(LIBIFPACK) \
           $(LIBLAPACK) $(LIBBLAS) $(LIBY12M) \
           $(LIBTRILINOS_UTIL)

#=======================================================================
# Epetra test source files
#=======================================================================

TEST_CC = vbr_main.cpp
TEST_C = 
TEST_F =

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC =

TEST_OBJ          =  $(TEST_CC:.cpp=.o) $(TEST_C:.c=.o)  $(TEST_F:.f=.o)

TARGET_MPI_MSR = vbr_msr_mpi
TARGET_SERIAL_MSR = vbr_msr_serial
TARGET_MPI_VBR = vbr_vbr_mpi
TARGET_SERIAL_VBR = vbr_vbr_serial
TARGET = $(TARGET_$(TRILINOS_COMM)_$(FORMAT))


$(TARGET): $(TEST_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.epetra

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI_MSR) $(TARGET_SERIAL_MSR) \
                    $(TARGET_MPI_VBR) $(TARGET_SERIAL_VBR)
