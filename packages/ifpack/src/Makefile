include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(IFPACK_ARCH_DEFINES) \
         $(EPETRA_ARCH_DEFINES) \
         $(AZTECOO_ARCH_DEFINES) -D$(EPETRA_COMM) -D$(AZTECOO_COMM) \
		 $(EPETRA_DEFAULT_DEFINES)

INCLUDES = $(ARCH_INCLUDES) $(IFPACK_INCLUDES) \
           $(EPETRA_INCLUDES) $(AZTECOO_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS = $(CXXFLAGS)
#=======================================================================
# IFPACK source files
#=======================================================================


IFPACK_CC_new = Ifpack_IlukGraph.cpp Ifpack_CrsRiluk.cpp
IFPACK_CC_old = \
            ifp_BlockMat.cc          ifp_BlockVec.cc \
            ifp_DenseMat.cc           \
            ifp_SparseUtil.cc        ifp_biluk.cc \
            ifp_brelax.cc            ifp_c_wrappers.cc \
            ifp_spmm.cc              ifp_spsm.cc

IFPACK_C  = az_ifpack_prec_create.c   az_ifpack_precon.c  \
            az_ifpack_iterate.c       az_ifpack_prec_destroy.c \
            az_ifpack_solve.c



#=======================================================================
# IFPACK include files
#=======================================================================

IFPACK_INC =

IFPACK_OBJ          = $(IFPACK_CC_new:.cpp=.o) $(IFPACK_CC_old:.cc=.o)  $(IFPACK_C:.c=.o) 

$(LIBIFPACK): $(IFPACK_OBJ)
	$(AR) $(ARFLAGS) $(LIBIFPACK) $(IFPACK_OBJ)
	$(RANLIB) $(LIBIFPACK)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.ifpack

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  $(LIBIFPACK) *~
