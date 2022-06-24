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
CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#set(MATLAB_FIND_DEBUG TRUE)

# --> ACTION: define Trilinos install directory
#set(TRILINOS_INSTALL_DIR "/Users/pakuber/Compadre/Trilinos-intrelab/pthread/install_release")
set(TRILINOS_INSTALL_DIR "/ascldap/users/pakuber/Compadre/Trilinos-intrelab/build/install")
if(DEFINED ENV{TRILINOS_INSTALL_DIR})
  set(TRILINOS_INSTALL_DIR $ENV{TRILINOS_INSTALL_DIR})
endif()


# --> ACTION: define Matlab root directory
#set(Matlab_ROOT_DIR "/Applications/MATLAB_R2020a.app/")
set(Matlab_ROOT_DIR "/usr/local/matlab/R2021b/")
if(DEFINED ENV{Matlab_ROOT_DIR})
  set(Matlab_ROOT_DIR $ENV{Matlab_ROOT_DIR})
endif()
find_package(Matlab REQUIRED)

# --> ACTION: define install directory (full path or relative to the intrelab directory)
set(INST_DIR "install")


set(CMAKE_PREFIX_PATH ${TRILINOS_INSTALL_DIR} ${CMAKE_PREFIX_PATH})
if(IS_ABSOLUTE ${INST_DIR})
  set(CMAKE_INSTALL_PREFIX ${INST_DIR}) 
else()
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/${INST_DIR}") 
endif()

FIND_PACKAGE(Trilinos REQUIRED)

set(MOBJS "intrepid_setJacobian;intrepid_setJacobianInv;intrepid_setJacobianDet;intrepid_getNumCubaturePoints;intrepid_getCubature;intrepid_getBasisValues;intrepid_computeCellMeasure;intrepid_HGRADtransformGRAD;intrepid_HGRADtransformVALUE;intrepid_multiplyMeasure;intrepid_integrate;intrepid_mapToPhysicalFrame;intrepid_mapToReferenceFrame;intrepid_scalarMultiplyDataField;intrepid_scalarMultiplyDataData;intrepid_dataIntegral;intrepid_evaluate;intrepid_multiplyVectorDataField;intrepid_dotMultiplyDataField;intrepid_vectorMultiplyDataField;intrepid_dotMultiplyDataData;intrepid_computeEdgeMeasure;intrepid_computeFaceMeasure;intrepid_getCellDeviatoricMat;intrepid_getCellVolumetricMat;intrepid_getPhysicalSideNormals;intrepid_mapToReferenceSubcell;intrepid_tensorMultiplyDataData;intrepid_tensorMultiplyDataField;defGrad;evaluateVectorField;evaluateVectorGradField;externalForceElastostaticsTL;geometricBmat;intrepid_det;intrepid_inverse;intrepid_transpose;materialBmat;strainEulerianLagrangian;stress;tangentStiffnessElastostaticsTL;tensorMultiplyDataField;myMatVecProductDataField;myMatVecProductDataData")

if (NOT APPLE)
  matlab_add_mex(NAME m2i_helpers
      SHARED
      SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/m2i_helpers.cpp
  )
  target_link_libraries(m2i_helpers Trilinos::all_libs)
endif()

function(make_mex fname)
  if (APPLE)
    matlab_add_mex(NAME ${fname}
        SHARED
        SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/m2i_helpers.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp
    )
  else()
    matlab_add_mex(NAME ${fname}
        SHARED
        SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp
    )
    target_link_libraries(${fname} m2i_helpers)
  endif()
  target_link_libraries(${fname} Trilinos::all_libs)
endfunction(make_mex)

foreach(fname ${MOBJS})
  make_mex(${fname})
endforeach()

install(TARGETS ${MOBJS}
  DESTINATION ${CMAKE_INSTALL_PREFIX} 
)
