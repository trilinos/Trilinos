include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM) \
         -DIFPACK -DSPARSEBLAS

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src \
           -I$(TRILINOS_HOME)/packages/aztecoo/src \
           -I$(TRILINOS_HOME)/packages/triutils/src

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)


FORMAT=MSR

LIB_PATHS= $(LIBPETRA) $(LIBAZTEC) $(LIBIFPACK) $(LIBMACHDEP) \
           $(LIBLAPACK) $(LIBSPBLAS) $(LIBBLAS) $(LIBY12M) \
           $(LIBTRIUTILS)

#=======================================================================
# Petra test source files
#=======================================================================

TEST_CC = cc_main.cc
TEST_C = 

TEST_F =

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = prototypes.h iohb.h

TEST_OBJ          = $(TEST_CC:.cc=.o) $(TEST_C:.c=.o) $(TEST_F:.f=.o)

TARGET_MPI_MSR = petra_transpose_mpi.exe
TARGET_SERIAL_MSR = petra_transpose_serial.exe
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
