# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBAZTECOO) by typing 'make'.
#

include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

FORMAT=MSR

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM) \
         -DIFPACK -DSPARSEBLAS -D$(FORMAT)

INCLUDES = $(ARCH_INCLUDES) $(EPETRA_INCLUDES) \
           $(AZTECOO_INCLUDES) $(TRIUTILS_INCLUDES) $(SPBLAS_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)



LIB_PATHS= $(LIBAZTECOO) $(LIBEPETRA) $(LIBIFPACK) \
           $(LIBLAPACK) $(LIBBLAS) $(LIBY12M) $(LIBSPBLAS) \
           $(LIBTRILINOS_UTIL)

#=======================================================================
# Epetra test source files
#=======================================================================

TEST_CC = 
TEST_C = \
          main.c \
          create_vbr.c \
          distrib_msr_matrix.c \
          distrib_vbr_matrix.c \
          iohb.c \
          read_coo.c \
          read_hb.c \
          scscmv.c \
          scscres.c \
          smsrres.c \
          svbrres.c \
          write_vec.c

TEST_F = blassm.f \
          formats.f \
          unary.f \

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC =

TEST_OBJ          =  $(TEST_CC:.cc=.o) $(TEST_C:.c=.o)  $(TEST_F:.f=.o)

TARGET_MPI_MSR = azoo_msr_mpi.exe
TARGET_SERIAL_MSR = azoo_msr_serial.exe
TARGET_MPI_VBR = azoo_vbr_mpi.exe
TARGET_SERIAL_VBR = azoo_vbr_serial.exe
TARGET = $(TARGET_$(TRILINOS_COMM)_$(FORMAT))


$(TARGET): $(TEST_OBJ)
	$(LINKER) $(LDFLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.epetra

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI_MSR) $(TARGET_SERIAL_MSR) \
                    $(TARGET_MPI_VBR) $(TARGET_SERIAL_VBR)
