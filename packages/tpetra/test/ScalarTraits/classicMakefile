# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#

TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

LIBPETRA= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libpetra.a

include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

# Petra communication defines
PETRA_COMM_SERIAL          = SERIAL
PETRA_COMM_MPI             = PETRA_MPI
PETRA_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))

TSF_COMM_SERIAL          = TSF_SERIAL
TSF_COMM_MPI             = TSF_MPI
TSF_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))


DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM) -D$(TSF_COMM)

INCLUDES = $(ARCH_INCLUDES) \
	   -I$(TRILINOS_HOME)/src/tpetra \
           -I$(TRILINOS_HOME)/src/petra \
           $(BLAS_INCLUDES) 

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
LDFLAGS=$(ARCH_LDFLAGS)


LIB_PATHS= $(LIBPETRA_TSF) $(LIBPETRA) $(LIBTPETRA) $(LIBLAPACK) $(LIBBLAS)

#=======================================================================
# Petra test source files
#=======================================================================

TEST_C =

TEST_CC = cc_main.cc

TEST_F = 
#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = 

TEST_C_OBJ          = $(TEST_C:.c=.o)
TEST_CC_OBJ          = $(TEST_CC:.cc=.o)
TEST_F_OBJ          = $(TEST_F:.f=.o)

TARGET_C_MPI = c_tpetra_densematrix_mpi
TARGET_CC_MPI = cc_tpetra_densematrix_mpi
TARGET_F_MPI = f_tpetra_densematrix_mpi
TARGET_C_SERIAL = c_tpetra_densematrix_serial
TARGET_CC_SERIAL = cc_tpetra_densematrix_serial
TARGET_F_SERIAL = f_tpetra_densematrix_serial

TARGET_C = $(TARGET_C_$(TRILINOS_COMM))
TARGET_CC = $(TARGET_CC_$(TRILINOS_COMM))
TARGET_F = $(TARGET_F_$(TRILINOS_COMM))

all: $(TARGET_C) $(TARGET_CC) $(TARGET_F)

c: $(TARGET_C)
cc: $(TARGET_CC)
f: $(TARGET_F)


$(TARGET_C): $(TEST_C_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_C_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI)  -o $(TARGET_C)

$(TARGET_CC): $(TEST_CC_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_CC_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI)  -o $(TARGET_CC)

$(TARGET_F): $(TEST_F_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_F_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI)  -o $(TARGET_F)

#
# dependencies for 'f' files (none at this time)
#
#include ../../etc/depends.petra

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_C_MPI) $(TARGET_C_SERIAL) \
                    $(TARGET_CC_MPI) $(TARGET_CC_SERIAL) \
                    $(TARGET_F_MPI) $(TARGET_F_SERIAL) 
                     
