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

LIBMACHDEP= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libmachdep.a
LIBIFPACK= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libifpack.a
LIBPETRA= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libpetra.a

include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

# Petra communication defines
PETRA_COMM_SERIAL          = SERIAL
PETRA_COMM_MPI             = PETRA_MPI
PETRA_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM) \
         -DIFPACK -DSPARSEBLAS
# -DMULTI_VECTOR

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/src/ifpack \
           -I$(TRILINOS_HOME)/src/petra $(BLAS_INCLUDES) \
           $(SPBLAS_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS) $(DEFINES)


FORMAT=MSR

LIB_PATHS= $(LIBIFPACK) $(LIBPETRA) $(LIBMACHDEP) $(LIBSPBLAS) $(LIBBLAS)

#=======================================================================
# Petra test source files
#=======================================================================

TEST_CC = cc_main.cc
TEST_C = read_hb.c write_vec.c distrib_msr_matrix.c smsrres.c \
           scscres.c read_coo.c scscmv.c distrib_vbr_matrix.c \
           create_vbr.c svbrres.c iohb.c

TEST_F = formats.f

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = prototypes.h iohb.h

TEST_OBJ          = $(TEST_CC:.cc=.o) $(TEST_C:.c=.o) $(TEST_F:.f=.o)

TARGET_MPI_MSR = CC_petra_hb_mpi_msr
TARGET_SERIAL_MSR = CC_petra_hb_serial_msr
TARGET_MPI_VBR = CC_petra_hb_mpi_vbr
TARGET_SERIAL_VBR = CC_petra_hb_serial_vbr
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
