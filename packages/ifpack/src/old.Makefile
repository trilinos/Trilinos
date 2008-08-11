# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBIFPACK) by typing 'make'.
#


TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

LIBIFPACK = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libifpack.a

include $(TRILINOS_HOME)/etc/makefile.$(TRILINOS_TARGET)

# IFPACK communication defines
IFPACK_COMM_SERIAL          = SERIAL
IFPACK_COMM_MPI             = AZTEC_MPI
IFPACK_COMM                 = $(IFPACK_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(IFPACK_ARCH_DEFINES) -D$(IFPACK_COMM) \
         -DIFPACK

INCLUDES = -I. $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/src/aztec $(BLAS_INCLUDES) 

LIB_PATHS=

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS = $(CXXFLAGS)
#=======================================================================
# IFPACK source files
#=======================================================================


IFPACK_CC = ifp_BlockMat.cc		  ifp_BlockVec.cc \
            ifp_DenseMat.cc		   \
            ifp_SparseUtil.cc		  ifp_biluk.cc \
            ifp_brelax.cc		  ifp_c_wrappers.cc \
            ifp_spmm.cc			  ifp_spsm.cc

IFPACK_C  = az_ifpack_prec_create.c   az_ifpack_precon.c  \
            az_ifpack_iterate.c       az_ifpack_prec_destroy.c \
            az_ifpack_solve.c


#=======================================================================
# IFPACK include files
#=======================================================================

IFPACK_INC = \
az_ifpack.h	   ifp_BlockMat.h      ifp_BlockVec.h \
ifp_DenseMat.h	   ifp_GlobalPrecon.h  ifp_LocalMat.h \
ifp_LocalPrecon.h  ifp_Matrix.h        ifp_Precon.h \
ifp_SparseUtil.h   ifp_arch.h	       ifp_biluk.h \
ifp_blas1.h	   ifp_blas3.h	       ifp_brelax.h \
ifp_c_wrappers.h   ifp_ifpack.h        ifp_lapackd.h \
ifp_spblas.h

IFPACK_OBJ          = $(IFPACK_CC:.cc=.o)  $(IFPACK_C:.c=.o)

$(LIBIFPACK): $(IFPACK_OBJ)
	$(AR) rcv $(LIBIFPACK) $(IFPACK_OBJ)
	$(RANLIB) $(LIBIFPACK)

#
# dependencies for 'f' files (none at this time)
#
#include ../../etc/depends.ifpack

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  $(LIBIFPACK) *~
