include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

# Epetra communication defines
EPETRA_MPICOMM_CLASS = Epetra_MpiComm.cpp Epetra_MpiDistributor.cpp
EPETRA_SERIALCOMM_CLASS = 
EPETRA_EXTRACOMM_CLASS = $(EPETRA_$(TRILINOS_COMM)COMM_CLASS)

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM) $(EPETRA_DEFAULT_DEFINES)

INCLUDES = $(ARCH_INCLUDES) $(EPETRA_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
#=======================================================================
# EPETRA source files
#=======================================================================


ALL_EPETRA_CC = \
            Epetra_BLAS.cpp \
            Epetra_BlockMap.cpp \
            Epetra_CrsGraph.cpp \
            Epetra_CrsMatrix.cpp \
            Epetra_DenseMatrix.cpp \
            Epetra_BasicDirectory.cpp \
            Epetra_EasySerialDenseProblem.cpp \
            Epetra_EasySerialSpdDenseProblem.cpp \
            Epetra_Export.cpp \
            Epetra_Flops.cpp \
            Epetra_HardSerialDenseProblem.cpp \
            Epetra_HardSerialSpdDenseProblem.cpp \
            Epetra_Import.cpp \
            Epetra_LAPACK.cpp \
            Epetra_LinearProblem.cpp \
            Epetra_LocalMap.cpp \
            Epetra_Map.cpp \
            Epetra_MpiComm.cpp \
            Epetra_MultiVector.cpp \
            Epetra_Object.cpp \
            Epetra_SerialComm.cpp \
            Epetra_SerialDenseMatrix.cpp \
            Epetra_SerialDenseVector.cpp \
            Epetra_SerialSpdDenseMatrix.cpp \
            Epetra_Time.cpp \
            Epetra_Util.cpp \
            Epetra_VbrMatrix.cpp \
            Epetra_Vector.cpp \
            Epetra_MpiDistributor.cpp \
            Epetra_SerialDistributor.cpp \
            Epetra_ReducedLinearProblem.cpp
#           Epetra_C_wrappers.cpp Epetra_Fortran_wrappers.cpp

EPETRA_CC = Epetra_SerialComm.cpp \
            Epetra_SerialDistributor.cpp \
            Epetra_BlockMap.cpp \
            Epetra_Map.cpp \
            Epetra_LocalMap.cpp \
            Epetra_BasicDirectory.cpp \
            Epetra_Export.cpp \
            Epetra_Import.cpp \
            Epetra_MultiVector.cpp \
            Epetra_Vector.cpp \
            Epetra_BLAS.cpp \
            Epetra_CrsGraph.cpp \
            Epetra_CrsMatrix.cpp \
            Epetra_VbrMatrix.cpp \
            Epetra_DistObject.cpp \
            Epetra_CompObject.cpp \
            Epetra_SerialDenseMatrix.cpp \
            Epetra_SerialDenseVector.cpp \
            Epetra_SerialDenseSolver.cpp\
            Epetra_SerialSymDenseMatrix.cpp \
            Epetra_SerialSpdDenseSolver.cpp\
            Epetra_LAPACK.cpp \
            Epetra_LinearProblem.cpp \
            Epetra_Time.cpp Epetra_Object.cpp Epetra_Util.cpp Epetra_Flops.cpp \
            Epetra_ReducedLinearProblem.cpp \
            $(EPETRA_EXTRACOMM_CLASS)

EPETRA_F =

#=======================================================================
# EPETRA include files
#=======================================================================

EPETRA_INC =

EPETRA_OBJ          = $(EPETRA_CC:.cpp=.o) $(EPETRA_C:.c=.o) $(EPETRA_F:.f=.o)

$(LIBEPETRA): $(EPETRA_OBJ)
	mkdir -p $(TRILINOS_LIB_PATH)
	$(AR) $(ARFLAGS) $(LIBEPETRA) $(EPETRA_OBJ) $(TEMPLATE_OBJS)
	$(RANLIB) $(LIBEPETRA)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.epetra

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  $(LIBEPETRA) *~
