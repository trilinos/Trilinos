# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBAZTEC) by typing 'make'.
#
include $(TRILINOS_HOME)/build/TARGET

TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

#LIBMACHDEP= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libmachdep.a
LIBEPETRA= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libepetra.a
LIBAZTEC= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libaztecoo.a
LIBIFPACK= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libifpack.a
LIBASCI_NLS= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libnlspack.a

include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

# Petra communication defines
PETRA_COMM_SERIAL          = SERIAL
PETRA_COMM_MPI             = EPETRA_MPI
PETRA_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM) \
         -DIFPACK -DWITHSPBLAS

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src $(BLAS_INCLUDES) \
           $(SPBLAS_INCLUDES) -I$(TRILINOS_HOME)/packages/aztecoo/src \
           -I$(TRILINOS_HOME)/packages/nlspack/ \
           -I$(TRILINOS_HOME)/packages/nlspack/test_linear 

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)


 LIB_PATHS= $(LIBASCI_NLS) $(LIBEPETRA) $(LIBAZTEC) $(LIBIFPACK) \
	    $(LIBMACHDEP) $(LIBLAPACK) $(LIBSPBLAS) $(LIBBLAS) $(LIBY12M) 

#=======================================================================
# Petra test source files
#=======================================================================

TEST_CC = cc_main.cc fill.cc basis.cc NLS_Interface.cc
TEST_C = 
TEST_F = 

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = 

TEST_OBJ          = $(TEST_CC:.cc=.o) $(TEST_C:.c=.o) $(TEST_F:.f=.o)


TARGET_MPI = test_mpi
TARGET_SERIAL = test_serial
TARGET = $(TARGET_$(TRILINOS_COMM))


$(TARGET): $(TEST_OBJ) $(LIBASCI_NLS)
	@echo "Making " $(TARGET)
	$(LINKER) $(ARCH_FLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI) $(TARGET_SERIAL)
