include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

FORMAT=NONE
VBRMATRIX_MSR =
VBRMATRIX_VBR = -DVBRMATRIX

DEFINES= -D$(TRILINOS_ARCH) $(AZTECOO_ARCH_DEFINES) -D$(AZTECOO_COMM) \
         -DIFPACK -DSPARSEBLAS  $(VBRMATRIX_$(FORMAT))

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/komplex/src \
           -I$(TRILINOS_HOME)/packages/ifpack/src \
           -I$(TRILINOS_HOME)/packages/aztecoo/src \
           -I$(TRILINOS_HOME)/external/spblas

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(CXXFLAGS)

LIB_PATHS= $(LIBKOMPLEX) $(LIBAZTECOO) $(LIBIFPACK) $(LIBY12M) $(LIBMACHDEP) \
           $(LIBLAPACK) $(LIBSPBLAS) $(LIBBLAS) $(LIBTRILINOS_UTIL)

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
	$(LINKER) $(ARCH_FLAGS) simple.o $(LIB_PATHS) $(ARCH_LIBS) $(LIBMPI) -o simple
	@echo " "
	@echo " "
	@echo "Simple test driver is built"
