
include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)



DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM)

INCLUDES = $(ARCH_INCLUDES) \
	   -I$(TRILINOS_HOME)/packages/tpetra/src \
           $(BLAS_INCLUDES) 

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
LDFLAGS=$(ARCH_LDFLAGS)


LIB_PATHS= $(LIBEPETRA) $(LIBLAPACK) $(LIBBLAS)

#=======================================================================
# Petra test source files
#=======================================================================

TEST_C =

TEST_CXX = cxx_main.cpp

TEST_F = 
#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = 

TEST_C_OBJ          = $(TEST_C:.c=.o)
TEST_CXX_OBJ          = $(TEST_CXX:.cpp=.o)
TEST_F_OBJ          = $(TEST_F:.f=.o)

TARGET_C_MPI = c_tpetra_commtest_mpi
TARGET_CXX_MPI = cxx_tpetra_commtest_mpi
TARGET_F_MPI = f_tpetra_commtest_mpi
TARGET_C_SERIAL = c_tpetra_commtest_serial
TARGET_CXX_SERIAL = cxx_tpetra_commtest_serial
TARGET_F_SERIAL = f_tpetra_commtest_serial

TARGET_C = $(TARGET_C_$(TRILINOS_COMM))
TARGET_CXX = $(TARGET_CXX_$(TRILINOS_COMM))
TARGET_F = $(TARGET_F_$(TRILINOS_COMM))

all: $(TARGET_CXX)

c: $(TARGET_C)
cc: $(TARGET_CXX)
f: $(TARGET_F)


$(TARGET_C): $(TEST_C_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_C_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI)  -o $(TARGET_C)

$(TARGET_CXX): $(TEST_CXX_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_CXX_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI)  -o $(TARGET_CXX)

$(TARGET_F): $(TEST_F_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_F_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI)  -o $(TARGET_F)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.petra

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_C_MPI) $(TARGET_C_SERIAL) \
                    $(TARGET_CXX_MPI) $(TARGET_CXX_SERIAL) \
                    $(TARGET_F_MPI) $(TARGET_F_SERIAL)

