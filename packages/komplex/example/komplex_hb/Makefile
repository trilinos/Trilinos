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

LIBKOMPLEX = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libkomplex.a
LIBAZTEC = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libaztec.a
LIBIFPACK= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libifpack.a
LIBMACHDEP= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libmachdep.a

include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

# Aztec communication defines
AZTEC_COMM_SERIAL          = SERIAL
AZTEC_COMM_MPI             = AZTEC_MPI
AZTEC_COMM                 = $(AZTEC_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(AZTEC_ARCH_DEFINES) -D$(AZTEC_COMM) \
         -DIFPACK -DSPARSEBLAS

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/src/komplex \
           -I$(TRILINOS_HOME)/src/ifpack \
           -I$(TRILINOS_HOME)/src/aztec $(BLAS_INCLUDES) \
           $(SPBLAS_INCLUDES)

FORMAT=NONE
VBRMATRIX_MSR =
VBRMATRIX_VBR = -DVBRMATRIX

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES) $(VBRMATRIX_$(FORMAT))
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(VBRMATRIX_$(FORMAT))
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES) $(VBRMATRIX_$(FORMAT))
LDFLAGS=$(CXXFLAGS)

LIB_PATHS= $(LIBKOMPLEX) $(LIBAZTEC) $(LIBIFPACK) $(LIBY12M) $(LIBMACHDEP) \
           $(LIBLAPACK) $(LIBSPBLAS) $(LIBBLAS)

#=======================================================================
# Komplex test source files
#=======================================================================

TEST_C = main.c read_hb.c distrib_msr_matrix.c smsrres.c \
           scscres.c read_coo.c scscmv.c distrib_vbr_matrix.c \
           create_vbr.c svbrres.c iohb.c

TEST_F = formats.f blassm.f unary.f

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = prototypes.h iohb.h

TEST_OBJ          = $(TEST_C:.c=.o) $(TEST_F:.f=.o)

TARGET_MPI_MSR = komplex_hb_mpi_msr
TARGET_SERIAL_MSR = komplex_hb_serial_msr
TARGET_MPI_VBR = komplex_hb_mpi_vbr
TARGET_SERIAL_VBR = komplex_hb_serial_vbr
TARGET = $(TARGET_$(TRILINOS_COMM)_$(FORMAT))


all: $(FORMAT)

# application code

$(TARGET): $(TEST_OBJ)
	$(LINKER) $(ARCH_FLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
     $(LIBMPI) -o $(TARGET)
clean:
	@echo "cleaning ..."
	/bin/rm -f $(TEST_OBJ)  $(TARGET_MPI_MSR) $(TARGET_SERIAL_MSR) \
                             $(TARGET_MPI_VBR) $(TARGET_SERIAL_VBR) \
                             simple.o simple

VBR: $(TARGET)
	@echo " "
	@echo " "
	@echo "Test Driver for VBR Matrices is built"

MSR: $(TARGET)
	@echo " "
	@echo " "
	@echo "Test Driver for MSR Matrices is built"

NONE:
	@echo " "
	@echo " "
	@echo "Please type one of the following:"
	@echo " "
	@echo "  make simple      (simple test driver)"
	@echo "  make FORMAT=MSR  (Harwell-Boeing test driver for MSR matrices)"
	@echo "  make FORMAT=VBR  (Harwell-Boeing test driver for VBR matrices)"
	@echo " "
	@echo "  NOTE:  You may need to use gmake if make gives errors"
	@echo " "

simple: simple.o
	$(LD_COMPILER) simple.o $(LDFLAGS) -o simple
	@echo " "
	@echo " "
	@echo "Simple test driver is built"
