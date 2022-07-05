#
# This CMakeLists.txt file configures the Intrepid / Matlab interface.
#
# This interface has dependencies on Trilinos packages.
#
# Note that the Trilinos install can be either a RELEASE
# or a DEBUG install; for performance, the RELEASE install
# is recommended.
#

PROJECT(Intrelab)
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(MATLAB_FIND_DEBUG TRUE)

# --> ACTIONS:
# 1.) Specify Trilinos install directory with -DTrilinos_ROOT=     (case-sensitive)
# 2.) Specify Matlab install directory with   -DMatlab_ROOT_DIR=   (case-sensitive and note extra _DIR)
# 3.) Specify directory to install intrelab   -DINST_DIR=          (case-sensitive)

# handle Matlab_ROOT_DIR and all variants 
# allow for {Matlab_ROOT_DIR, MATLAB_ROOT_DIR, Matlab_ROOT, MATLAB_ROOT}
if (MATLAB_ROOT_DIR)
  if (NOT Matlab_ROOT_DIR)
    set(Matlab_ROOT_DIR ${MATLAB_ROOT_DIR})
  endif()
endif()
if (Matlab_ROOT)
  if (NOT Matlab_ROOT_DIR)
    set(Matlab_ROOT_DIR ${Matlab_ROOT})
  endif()
endif()
if (MATLAB_ROOT)
  if (NOT Matlab_ROOT_DIR)
    set(Matlab_ROOT_DIR ${MATLAB_ROOT})
  endif()
endif()
# look for environment variable and give it precedence
# allow for {MATLAB_ROOT, Matlab_ROOT, MATLAB_ROOT_DIR, Matlab_ROOT_DIR}
if(DEFINED ENV{MATLAB_ROOT})
  set(Matlab_ROOT_DIR $ENV{MATLAB_ROOT})
endif()
if(DEFINED ENV{Matlab_ROOT})
  set(Matlab_ROOT_DIR $ENV{Matlab_ROOT})
endif()
if(DEFINED ENV{MATLAB_ROOT_DIR})
  set(Matlab_ROOT_DIR $ENV{MATLAB_ROOT_DIR})
endif()
if(DEFINED ENV{Matlab_ROOT_DIR})
  set(Matlab_ROOT_DIR $ENV{Matlab_ROOT_DIR})
endif()
if (NOT Matlab_ROOT_DIR)
  message(WARNING "Matlab_ROOT_DIR has NOT been specified.")    
endif()

find_package(Matlab REQUIRED)

# handle Trilinos_ROOT and all variants 
# allow for {Trilinos_ROOT, TRILINOS_ROOT}
if (TRILINOS_ROOT)
  if (NOT Trilinos_ROOT)
    set(Trilinos_ROOT ${TRILINOS_ROOT})
  endif()
endif()
# look for environment variable and give it precedence
# allow for {Trilinos_ROOT, TRILINOS_ROOT}
if(DEFINED ENV{TRILINOS_ROOT})
  set(Trilinos_ROOT $ENV{TRILINOS_ROOT})
endif()
if(DEFINED ENV{Trilinos_ROOT})
  set(Trilinos_ROOT $ENV{Trilinos_ROOT})
endif()
if (NOT Trilinos_ROOT)
  message(WARNING "Trilinos_ROOT has NOT been specified.")    
endif()

FIND_PACKAGE(Trilinos REQUIRED)

# --> ACTION: define install directory (full path or relative to the intrelab directory)
set(INST_DIR "install")

set(CMAKE_PREFIX_PATH ${TRILINOS_INSTALL_DIR} ${CMAKE_PREFIX_PATH})
if(IS_ABSOLUTE ${INST_DIR})
  set(CMAKE_INSTALL_PREFIX ${INST_DIR}) 
else()
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/${INST_DIR}") 
endif()

set(MOBJS intrepid_setJacobian intrepid_setJacobianInv intrepid_setJacobianDet
    intrepid_getNumCubaturePoints intrepid_getCubature intrepid_getBasisValues
    intrepid_computeCellMeasure intrepid_HGRADtransformGRAD
    intrepid_HGRADtransformVALUE intrepid_multiplyMeasure intrepid_integrate
    intrepid_mapToPhysicalFrame intrepid_mapToReferenceFrame
    intrepid_scalarMultiplyDataField intrepid_scalarMultiplyDataData
    intrepid_dataIntegral intrepid_evaluate intrepid_multiplyVectorDataField
    intrepid_dotMultiplyDataField intrepid_vectorMultiplyDataField
    intrepid_dotMultiplyDataData intrepid_computeEdgeMeasure
    intrepid_computeFaceMeasure intrepid_getCellDeviatoricMat
    intrepid_getCellVolumetricMat intrepid_getPhysicalSideNormals
    intrepid_mapToReferenceSubcell intrepid_tensorMultiplyDataData
    intrepid_tensorMultiplyDataField defGrad;evaluateVectorField
    evaluateVectorGradField externalForceElastostaticsTL geometricBmat
    intrepid_det intrepid_inverse intrepid_transpose materialBmat
    strainEulerianLagrangian stress tangentStiffnessElastostaticsTL
    tensorMultiplyDataField myMatVecProductDataField myMatVecProductDataData)

function(make_mex fname)
  matlab_add_mex(NAME ${fname}
      SHARED
      SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/m2i_helpers.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp
  )
  target_link_libraries(${fname} Trilinos::all_libs)
endfunction()

foreach(fname ${MOBJS})
  make_mex(${fname})
endforeach()

install(TARGETS ${MOBJS}
  DESTINATION ${CMAKE_INSTALL_PREFIX} 
)
