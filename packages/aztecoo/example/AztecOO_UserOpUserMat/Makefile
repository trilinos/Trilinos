include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM) \
         -DIFPACK -DSPARSEBLAS

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src \
           -I$(TRILINOS_HOME)/packages/aztecoo/src \
           -I$(TRILINOS_HOME)/packages/triutils/src

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)



LIB_PATHS= $(LIBAZTECOO) $(LIBEPETRA) $(LIBIFPACK) \
           $(LIBLAPACK) $(LIBBLAS) $(LIBY12M) \
           $(LIBTRIUTILS)

#=======================================================================
# Epetra test source files
#=======================================================================

TEST_CC = cxx_main.cpp Poisson2dOperator.cpp
TEST_C = 
TEST_F =

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC =

TEST_OBJ          =  $(TEST_CC:.cpp=.o) $(TEST_C:.c=.o)  $(TEST_F:.f=.o)

TARGET_MPI = poisson2d_mpi.exe
TARGET_SERIAL = poisson2d_serial.exe
TARGET = $(TARGET_$(TRILINOS_COMM))


$(TARGET): $(TEST_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI) $(TARGET_SERIA)
