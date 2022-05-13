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

# --> ACTION: define Trilinos install directory
set(TRILINOS_INSTALL_DIR "/Users/pakuber/Compadre/Trilinos-new/pthread/install_release")
set(CMAKE_PREFIX_PATH ${TRILINOS_INSTALL_DIR} ${CMAKE_PREFIX_PATH})

# --> ACTION: define Matlab binary directory
set(MATLAB_BIN_DIR "/Applications/MATLAB_R2020a.app/bin/")

# --> ACTION: define install directory (full path or relative to the intrelab directory)
set(INST_DIR "install")
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/${INST_DIR}") 


FIND_PACKAGE(Trilinos REQUIRED)

# set Matlab mex path with -largeArrayDims option
set(MEX "${MATLAB_BIN_DIR}/mex")
message(STATUS "MEX: ${MEX}")

# set compilers to Trilinos compilers
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})

# set flags to Trilinos flags
set(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})

# Trilinos doesn't propagate its C++ standard
set(CMAKE_CXX_STANDARD 17)
set(MEX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
if (APPLE)
  # Trilinos doesn't propagate its sysroot
  message(STATUS "${CMAKE_OSX_SYSROOT} Minimum OS X deployment version: ${CMAKE_OSX_DEPLOYMENT_TARGET}")
  set(MEX_FLAGS "--sysroot=${CMAKE_OSX_SYSROOT} -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET} ${MEX_FLAGS}")
endif(APPLE)


if (Trilinos_INCLUDE_DIRS) 
    set(Trilinos_INCLUDE_DIRS "-I${Trilinos_INCLUDE_DIRS}")
endif()
if (Trilinos_TPL_INCLUDE_DIRS) 
    set(Trilinos_TPL_INCLUDE_DIRS "-I${Trilinos_TPL_INCLUDE_DIRS}")
endif()
if (Trilinos_LIBRARY_DIRS) 
    set(Trilinos_LIBRARY_DIRS "-L${Trilinos_LIBRARY_DIRS}")
endif()
if (Trilinos_TPL_LIBRARY_DIRS) 
    set(Trilinos_TPL_LIBRARY_DIRS "-L${Trilinos_TPL_LIBRARY_DIRS}")
endif()

set(l_Trilinos "")
foreach(alib ${Trilinos_LIBRARIES})
  set(l_Trilinos "${l_Trilinos} -l${alib}")
endforeach()

set(l_Trilinos_TPL "")
foreach(alib ${Trilinos_TPL_LIBRARIES})
  if (APPLE)
    string(FIND ${alib} ".tbd" alib_result)
    # make sure .tbd not found in library name
  else()
    set(alib_result -1)
  endif(APPLE)
  if (alib_result LESS 0)
    set(l_Trilinos_TPL "${l_Trilinos_TPL} -l${alib}")
  else() 
    get_filename_component(libname ${alib} NAME_WLE)
    string(SUBSTRING ${libname} 3 -1 reduced_libname)
    set(l_Trilinos_TPL "${l_Trilinos_TPL} -l${reduced_libname}")
  endif()
endforeach()

separate_arguments(ll UNIX_COMMAND ${l_Trilinos})
set(l_Trilinos ${ll})

separate_arguments(ll UNIX_COMMAND ${l_Trilinos_TPL})
set(l_Trilinos_TPL ${ll})

if(APPLE)
  set(MEX_EXTENSION "mexmaci64")
else()
  set(MEX_EXTENSION "mexa64")
endif(APPLE)

# create objs folder
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/objs)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o
    COMMAND "${MEX}" "CXX=${CMAKE_CXX_COMPILER}" "CXXFLAGS=${MEX_FLAGS}" "-largeArrayDims" -c ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/src/m2i_helpers.cpp"
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/m2i_helpers.cpp
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/objs)

function(make_mex fname)

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o
      COMMAND "${MEX}" "CXX=${CMAKE_CXX_COMPILER}" "CXXFLAGS=${MEX_FLAGS}" "-largeArrayDims" -c ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp"
      BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/objs
      )
  
  add_custom_command(OUTPUT ${fname}.${MEX_EXTENSION}
      COMMAND "${MEX}" "CXX=${CMAKE_CXX_COMPILER}" "CXXFLAGS=${MEX_FLAGS}" "-largeArrayDims" ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} ${l_Trilinos} ${l_Trilinos_TPL} "${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o" "${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o"
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o ${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o
      BYPRODUCTS ${fname}.${MEX_EXTENSION}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
  
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${fname}.${MEX_EXTENSION} DESTINATION ${CMAKE_INSTALL_PREFIX})

endfunction(make_mex)

set(MOBJS "intrepid_setJacobian;intrepid_setJacobianInv;intrepid_setJacobianDet;intrepid_getNumCubaturePoints;intrepid_getCubature;intrepid_getBasisValues;intrepid_computeCellMeasure;intrepid_HGRADtransformGRAD;intrepid_HGRADtransformVALUE;intrepid_multiplyMeasure;intrepid_integrate;intrepid_mapToPhysicalFrame;intrepid_mapToReferenceFrame;intrepid_scalarMultiplyDataField;intrepid_scalarMultiplyDataData;intrepid_dataIntegral;intrepid_evaluate;intrepid_multiplyVectorDataField;intrepid_dotMultiplyDataField;intrepid_vectorMultiplyDataField;intrepid_dotMultiplyDataData;intrepid_computeEdgeMeasure;intrepid_computeFaceMeasure;intrepid_getCellDeviatoricMat;intrepid_getCellVolumetricMat;intrepid_getPhysicalSideNormals;intrepid_mapToReferenceSubcell;intrepid_tensorMultiplyDataData;intrepid_tensorMultiplyDataField;defGrad;evaluateVectorField;evaluateVectorGradField;externalForceElastostaticsTL;geometricBmat;intrepid_det;intrepid_inverse;intrepid_transpose;materialBmat;strainEulerianLagrangian;stress;tangentStiffnessElastostaticsTL;tensorMultiplyDataField;myMatVecProductDataField;myMatVecProductDataData")

# dummy targets to build up dependencies
add_custom_target(m2i_helpers
    ALL
    COMMAND # left intentionally blank
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o
    )
foreach(fname ${MOBJS})
    make_mex(${fname})
    # builds dependency chain
    add_custom_target(${fname}
        ALL
        COMMAND # left intentionally blank
        #DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${fname}.${MEX_EXTENSION}
        DEPENDS m2i_helpers ${fname}.${MEX_EXTENSION}
        )
endforeach()

if(APPLE)
  message("\"*.dylib was built for newer macOS version (*.*) than being linked (10.9)\" type warnings can be safely ignored. MEX compiler specifies a -mmacosx-min-version=10.9, which can not be overridden.")
endif(APPLE)
