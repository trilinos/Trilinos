include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

LIBSPARSE=$(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libepetra_sparse.a
SPARSE_INCLUDES=-I$(TRILINOS_HOME)/contrib/epetra_sparse/src


DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM)
INCLUDES = $(ARCH_INCLUDES) $(EPETRA_INCLUDES) \
           $(SPARSE_INCLUDES) \
           $(TRIUTILS_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)


FORMAT=MSR

LIB_PATHS= $(LIBEPETRA) $(LIBSPARSE) \
           $(LIBLAPACK) $(LIBBLAS) $(LIBY12M) \
           $(LIBTRIUTILS)

#=======================================================================
# Epetra test source files
#=======================================================================

TEST_CC = cpp_main.cpp
TEST_C = 

TEST_F =

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC =

TEST_OBJ          = $(TEST_CC:.cpp=.o) $(TEST_C:.c=.o) $(TEST_F:.f=.o)

TARGET_MPI_MSR = epetra_sparse_mpi
TARGET_SERIAL_MSR = epetra_sparse_serial
TARGET = $(TARGET_$(TRILINOS_COMM)_$(FORMAT))


$(TARGET): $(TEST_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI_MSR) $(TARGET_SERIAL_MSR) \
                    $(TARGET_MPI_VBR) $(TARGET_SERIAL_VBR)
