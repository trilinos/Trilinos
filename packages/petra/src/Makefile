include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM)

INCLUDES = $(ARCH_INCLUDES) $(PETRA_INCLUDES) $(AZTECOO_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
#=======================================================================
# PETRA source files
#=======================================================================


PETRA_CC = Aztec_OO.cc \
           Petra_Object.cc Petra_Comm.cc Petra_Time.cc \
           Petra_BlockMap.cc \
           Petra_Map.cc Petra_LocalMap.cc \
           Petra_RDP_MultiVector.cc  \
            Petra_RDP_Vector.cc \
            Petra_RDP_CRS_Matrix.cc \
            Petra_CRS_Graph.cc \
            Petra_Directory.cc \
            Petra_Import.cc \
            Petra_Export.cc \
            Petra_Flops.cc \
            Petra_Util.cc \
            Petra_BLAS.cc \
            Petra_LAPACK.cc \
            Petra_RDP_DenseMatrix.cc \
            Petra_RDP_DenseVector.cc \
            GSComm_Plan.cc \
            GSComm_Comm.cc \
            Petra_RDP_VBR_Matrix.cc \
            Petra_RDP_LinearProblem.cc
#           Petra_RDP_SPD_DenseMatrix.cc \
#           Petra_c_wrappers.cc Petra_fortran_wrappers.cc

PETRA_F =

#=======================================================================
# PETRA include files
#=======================================================================

PETRA_INC =

PETRA_OBJ = $(PETRA_CC:.cc=.$(OBJ_EXT)) $(PETRA_C:.c=.$(OBJ_EXT)) $(PETRA_F:.f=.$(OBJ_EXT))

#$(LIBPETRA): $(PETRA_OBJ)
LIBPETRA: $(PETRA_OBJ)
	mkdir -p $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)
	$(AR) $(ARFLAGS) $(LIB_WIN_FLAGS)$(LIBPETRA) $(PETRA_OBJ)
# following ranlib commented out by rbl
#	$(RANLIB) $(LIBPETRA)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.petra

clean:
	@echo "junk" > dummy.$(OBJ_EXT)
	@rm -f *.$(OBJ_EXT)  $(LIBPETRA) *~
