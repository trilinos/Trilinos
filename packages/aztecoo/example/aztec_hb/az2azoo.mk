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
LIBPETRA= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libpetra.a
LIBAZTEC= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libaztec.a
LIBIFPACK= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libifpack.a
LIBTRILINOS_UTIL= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libtriutils.a

include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

# Petra communication defines
PETRA_COMM_SERIAL          = SERIAL
PETRA_COMM_MPI             = PETRA_MPI
PETRA_COMM                 = $(PETRA_COMM_$(TRILINOS_COMM))

FORMAT=VBR

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM) \
         -DIFPACK -D$(FORMAT)

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/src/petra $(BLAS_INCLUDES) \
           -I$(TRILINOS_HOME)/src/aztec \
           -I$(TRILINOS_HOME)/src/triutils

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)



LIB_PATHS= $(LIBPETRA) $(LIBAZTEC) $(LIBIFPACK) $(LIBMACHDEP) \
           $(LIBLAPACK) $(LIBBLAS) $(LIBY12M) \
           $(LIBTRILINOS_UTIL)

#=======================================================================
# Petra test source files
#=======================================================================

TEST_C = \
az2azoo.c \
blassm.f \
create_vbr.c \
distrib_msr_matrix.c \
distrib_vbr_matrix.c \
formats.f \
iohb.c \
read_coo.c \
read_hb.c \
scscmv.c \
scscres.c \
smsrres.c \
svbrres.c \
unary.f \
write_vec.c

TEST_F =

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC =

TEST_OBJ          =  $(TEST_CC:.cc=.o) $(TEST_C:.c=.o)  $(TEST_F:.f=.o)

TARGET_MPI_MSR = az2azoo_oo_msr_mpi
TARGET_SERIAL_MSR = az2azoo_oo_msr_serial
TARGET_MPI_VBR = az2azoo_oo_vbr_mpi
TARGET_SERIAL_VBR = az2azoo_oo_vbr_serial
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
