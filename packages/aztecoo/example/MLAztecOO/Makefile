include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

FORMAT=MSR

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM) \
         -DIFPACK -DSPARSEBLAS -D$(FORMAT) -DML_MPI -DAZTEC

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src \
           -I$(TRILINOS_HOME)/packages/aztecoo/src \
           -I$(TRILINOS_HOME)/packages/triutils/src \
           $(EPETRA_INCLUDES) $(AZTECOO_INCLUDES) $(TRIUTILS_INCLUDES)  \
           -I/home/jhu/ml/trilinos/ml/Obj 

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)



LIB_PATHS= $(LIBAZTECOO) $(LIBIFPACK) $(LIBEPETRA) \
           $(LIBY12M) $(LIBTRIUTILS) \
           $(LIBML) $(LIBSUPERLU) \
           $(LIBLAPACK) $(LIBBLAS)

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

TEST_OBJ =  $(TEST_CC:.cpp=.$(OBJ_EXT)) $(TEST_C:.c=.$(OBJ_EXT))  $(TEST_F:.f=.$(OBJ_EXT))

TARGET_MPI_MSR = aztecoo_msr_mpi.exe
TARGET_SERIAL_MSR = aztecoo_msr_serial.exe
TARGET_MPI_VBR = aztecoo_vbr_mpi.exe
TARGET_SERIAL_VBR = aztecoo_vbr_serial.exe
TARGET = $(TARGET_$(TRILINOS_COMM)_$(FORMAT))


$(TARGET): $(TEST_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.epetra

clean:
	@echo "junk" > dummy.$(OBJ_EXT)
	@rm -f *.$(OBJ_EXT)  *~ $(TARGET_MPI_MSR) $(TARGET_SERIAL_MSR) \
                    $(TARGET_MPI_VBR) $(TARGET_SERIAL_VBR)
