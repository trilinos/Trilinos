# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBAZTEC) by typing 'make'.
#


TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

LIBPETRA= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libpetra.a
LIBTRILINOS_UTIL= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libtriutils.a
LIBSUPERLU=/u/mheroux/SuperLU/superlu_linux.a \
           /u/mheroux/SuperLU/blas_linux.a

include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

# Petra communication defines
PETRA_COMM_SERIAL          = SERIAL
PETRA_COMM_MPI             = PETRA_MPI
PETRA_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM)

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/src/petra \
           -I/u/mheroux/SuperLU/SRC \
           -I$(TRILINOS_HOME)/src/triutils

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)


FORMAT=MSR

LIB_PATHS= $(LIBPETRA) $(LIBSUPERLU) \
           $(LIBLAPACK) $(LIBBLAS) $(LIBY12M) \
           $(LIBTRILINOS_UTIL)

#=======================================================================
# Petra test source files
#=======================================================================

TEST_CC = cpp_main.cpp PetraSuperLU.cpp
TEST_C = 

TEST_F =

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = prototypes.h iohb.h

TEST_OBJ          = $(TEST_CC:.cpp=.o) $(TEST_C:.c=.o) $(TEST_F:.f=.o)

TARGET_MPI_MSR = petrasuperlu_mpi
TARGET_SERIAL_MSR = petrasuperlu_serial
TARGET = $(TARGET_$(TRILINOS_COMM)_$(FORMAT))


$(TARGET): $(TEST_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

#
# dependencies for 'f' files (none at this time)
#
#include ../../etc/depends.petra

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI_MSR) $(TARGET_SERIAL_MSR) \
                    $(TARGET_MPI_VBR) $(TARGET_SERIAL_VBR)
