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

LIBAZTEC= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libaztec.a
LIBIFPACK= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libifpack.a
LIBMACHDEP= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libmachdep.a

include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

# Petra communication defines
PETRA_COMM             = PETRA_MPI

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM) \
         -DAZ_MPI -DAZTEC_MPI -DIFPACK -DSPARSEBLAS

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/src/aztec $(BLAS_INCLUDES) 

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
LDFLAGS=$(ARCH_LDFLAGS)


LIB_PATHS= $(LIBAZTEC) $(LIBIFPACK) $(LIBY12M) $(LIBMACHDEP) $(LIBLAPACK) $(LIBSPBLAS) $(LIBBLAS) 

#=======================================================================
# Petra test source files
#=======================================================================

TEST_C = az_sync.c

TEST_CC = cc_main.cc

TEST_F = f_main.f
#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = 

TEST_CC_OBJ          = $(TEST_CC:.cc=.o)
TEST_C_OBJ          = $(TEST_C:.c=.o)

TARGET = az_sync_test

all: $(TARGET)
$(TARGET): $(TEST_CC_OBJ) $(TEST_C_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_CC_OBJ)  $(TEST_C_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI)  -o $(TARGET)

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET)