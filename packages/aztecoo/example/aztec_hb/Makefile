include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

# AztecOO communication defines
AZTEC_COMM_SERIAL          = SERIAL
AZTEC_COMM_MPI             = AZTEC_MPI
AZTEC_COMM                 = $(AZTEC_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(AZTEC_ARCH_DEFINES) -D$(AZTEC_COMM) \
         -DIFPACK -DSPARSEBLAS -DTEST_SINGULAR
#         -DSPARSEBLAS -DAZ_FUNKY_READ_UPDATE

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/ifpack/src \
           -I$(TRILINOS_HOME)/packages/aztecoo/src $(BLAS_INCLUDES) \
           $(SPBLAS_INCLUDES)

FORMAT=MSR
VBRMATRIX_MSR =
VBRMATRIX_VBR = -DVBRMATRIX

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES) $(VBRMATRIX_$(FORMAT))
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(VBRMATRIX_$(FORMAT))
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES) $(VBRMATRIX_$(FORMAT))
LDFLAGS=$(ARCH_LDFLAGS)


LIB_PATHS= $(LIBAZTECOO) $(LIBIFPACK) $(LIBY12M) $(LIBMACHDEP) \
           $(LIBLAPACK) $(LIBSPBLAS) $(LIBBLAS)

#=======================================================================
# AztecOO test source files
#=======================================================================

TEST_C = main.c read_hb.c write_vec.c distrib_msr_matrix.c smsrres.c \
           scscres.c read_coo.c scscmv.c distrib_vbr_matrix.c \
           create_vbr.c svbrres.c iohb.c
# az_read_update_funky.c

TEST_F = formats.f blassm.f unary.f

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = prototypes.h iohb.h

TEST_OBJ          = $(TEST_C:.c=.o) $(TEST_F:.f=.o)


TARGET_MPI_MSR = aztec_hb_mpi_msr.exe
TARGET_SERIAL_MSR = aztec_hb_serial_msr.exe
TARGET_MPI_VBR = aztec_hb_mpi_vbr.exe
TARGET_SERIAL_VBR = aztec_hb_serial_vbr.exe
TARGET = $(TARGET_$(TRILINOS_COMM)_$(FORMAT))


$(TARGET): $(TEST_OBJ)
	$(LINKER) $(LDFLAGS) $(ARCH_FLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.aztecoo

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI_MSR) $(TARGET_SERIAL_MSR) \
                    $(TARGET_MPI_VBR) $(TARGET_SERIAL_VBR)
