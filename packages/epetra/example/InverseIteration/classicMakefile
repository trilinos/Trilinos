include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

FORMAT=MSR

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM) \
         -DIFPACK -DSPARSEBLAS -D$(FORMAT)

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src \
           -I$(TRILINOS_HOME)/packages/ifpack/src \
           -I$(TRILINOS_HOME)/packages/triutils/src

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)



LIB_PATHS= $(LIBAZTECOO) $(LIBEPETRA) $(LIBIFPACK) \
           $(LIBLAPACK) $(LIBBLAS) $(LIBY12M) $(LIBSPBLAS) \
           $(LIBTRIUTILS)

#=======================================================================
# Epetra test source files
#=======================================================================

TEST_CC = cxx_main.cpp
TEST_C = 
TEST_F =

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC =

TEST_OBJ          =  $(TEST_CC:.cpp=.o) $(TEST_C:.c=.o)  $(TEST_F:.f=.o)

TARGET_MPI_MSR = inviter_msr_mpi.exe
TARGET_SERIAL_MSR = inviter_msr_serial.exe
TARGET_MPI_VBR = inviter_vbr_mpi.exe
TARGET_SERIAL_VBR = inviter_vbr_serial.exe
TARGET = $(TARGET_$(TRILINOS_COMM)_$(FORMAT))


$(TARGET): $(TEST_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

OPT_CXXFLAGS = -O3 #Put any other compiler flags here that you want

myfunction.o: myfunction.cpp
	$(CXX) $(OPT_CXXFLAGS) $(DEFINES) $(INCLUDES) -c myfunction.cpp #note that a tab must by first

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.epetra

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI_MSR) $(TARGET_SERIAL_MSR) \
                    $(TARGET_MPI_VBR) $(TARGET_SERIAL_VBR)
