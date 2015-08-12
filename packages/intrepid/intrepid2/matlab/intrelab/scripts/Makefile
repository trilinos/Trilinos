#
# This makefile compiles the Intrepid / Matlab interface.
#
# This interface has dependencies on Trilinos packages.
#
# The file Makefile.export.Trilinos_install is used to
# pull in all Trilinos dependencies.
#
# Point to your own Makefile.export.Trilinos_install,
# and point to your own MEX compiler.
#
# Note that the Trilinos install can be either a RELEASE
# or a DEBUG install; for performance, the RELEASE install
# is recommended.
#

# --> ACTION: include makefile export files ('install mode') for Trilinos; use RELEASE install
include /home/dridzal/research/Trilinos/BUILD/gcc-4.5.1_intrepid-ML_release/Makefile.export.Trilinos_install

# --> ACTION: point to your mex compiler; use the mex option -largeArrayDims 
MEX = /usr/local/MATLAB/R2010b/bin/mex -largeArrayDims

# --> ACTION: define install directory (full path or relative to the intrelab directory)
INST_DIR = install

# Set compilers to Trilinos compilers
CXX = $(Trilinos_CXX_COMPILER)
CPP = $(Trilinos_C_COMPILER)

# Set flags to Trilinos flags
CXXFLAGS = "$(Trilinos_CXX_COMPILER_FLAGS)"

# List target objects for compilation
MOBJS = intrepid_setJacobian.mexa64 \
        intrepid_setJacobianInv.mexa64 \
        intrepid_setJacobianDet.mexa64 \
        intrepid_getNumCubaturePoints.mexa64 \
        intrepid_getCubature.mexa64 \
        intrepid_getBasisValues.mexa64 \
        intrepid_computeCellMeasure.mexa64 \
        intrepid_HGRADtransformGRAD.mexa64 \
        intrepid_HGRADtransformVALUE.mexa64 \
        intrepid_multiplyMeasure.mexa64 \
        intrepid_integrate.mexa64 \
        intrepid_mapToPhysicalFrame.mexa64 \
        intrepid_mapToReferenceFrame.mexa64 \
        intrepid_scalarMultiplyDataField.mexa64 \
        intrepid_scalarMultiplyDataData.mexa64 \
        intrepid_dataIntegral.mexa64 \
        intrepid_evaluate.mexa64 \
        intrepid_multiplyVectorDataField.mexa64 \
        intrepid_dotMultiplyDataField.mexa64 \
        intrepid_vectorMultiplyDataField.mexa64 \
        intrepid_dotMultiplyDataData.mexa64 \
       	intrepid_computeEdgeMeasure.mexa64 \
        intrepid_computeFaceMeasure.mexa64 \
       	intrepid_getCellDeviatoricMat.mexa64 \
       	intrepid_getCellVolumetricMat.mexa64 \
       	intrepid_getPhysicalSideNormals.mexa64 \
       	intrepid_mapToReferenceSubcell.mexa64 \
       	intrepid_tensorMultiplyDataData.mexa64 \
       	intrepid_tensorMultiplyDataField.mexa64 \
	defGrad.mexa64 \
	evaluateVectorField.mexa64 \
	evaluateVectorGradField.mexa64 \
	externalForceElastostaticsTL.mexa64 \
	geometricBmat.mexa64 \
	intrepid_det.mexa64 \
	intrepid_inverse.mexa64 \
	intrepid_transpose.mexa64 \
	materialBmat.mexa64 \
	strainEulerianLagrangian.mexa64 \
	stress.mexa64 \
	tangentStiffnessElastostaticsTL.mexa64 \
	tensorMultiplyDataField.mexa64 \
        mlmex.mexa64 \
        myMatVecProductDataField.mexa64 \
        myMatVecProductDataData.mexa64

default: mex

mex: $(MOBJS)

# Use this if object files are needed:
#.SECONDARY:

%.o : src/%.cpp
	$(MEX) \
        CXX=$(CXX) LD=$(CXX) \
        CXXDEBUGFLAGS=$(CXXFLAGS) CXXOPTIMFLAGS=$(CXXFLAGS) -c \
        $(Trilinos_INCLUDE_DIRS) \
        src/$*.cpp

%.mexa64 : %.o m2i_helpers.o
	$(MEX) \
        CXX=$(CXX) LD=$(CXX) \
        CXXDEBUGFLAGS=$(CXXFLAGS) CXXOPTIMFLAGS=$(CXXFLAGS) LDDEBUGFLAGS=$(CXXFLAGS) LDOPTIMFLAGS=$(CXXFLAGS)  -cxx \
        $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) \
        $*.o m2i_helpers.o

install: mex
	mkdir -p $(INST_DIR); \
        mv *.mexa64 $(INST_DIR); \
        cp src/*.m $(INST_DIR)

clean:
	rm -f *.o *.mexa64; \
        rm -f $(INST_DIR)/*; \
        mkdir -p $(INST_DIR); \
        rmdir $(INST_DIR)

