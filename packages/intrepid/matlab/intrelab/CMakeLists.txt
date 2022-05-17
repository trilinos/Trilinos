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
#set(TRILINOS_INSTALL_DIR "/Users/pakuber/Compadre/Trilinos-new/pthread/install_release")
set(TRILINOS_INSTALL_DIR "/ascldap/users/pakuber/Compadre/Trilinos-intrelab/build/install")
if(DEFINED ENV{TRILINOS_INSTALL_DIR})
  set(TRILINOS_INSTALL_DIR $ENV{TRILINOS_INSTALL_DIR})
endif()


# --> ACTION: define Matlab binary directory
#set(MATLAB_BIN_DIR "/Applications/MATLAB_R2020a.app/bin/")
set(MATLAB_BIN_DIR "/usr/local/matlab/R2021b/bin/")
if(DEFINED ENV{MATLAB_BIN_DIR})
  set(MATLAB_BIN_DIR $ENV{MATLAB_BIN_DIR})
endif()

# --> ACTION: define install directory (full path or relative to the intrelab directory)
set(INST_DIR "install")


set(CMAKE_PREFIX_PATH ${TRILINOS_INSTALL_DIR} ${CMAKE_PREFIX_PATH})
if(IS_ABSOLUTE ${INST_DIR})
  set(CMAKE_INSTALL_PREFIX ${INST_DIR}) 
else()
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/${INST_DIR}") 
endif()

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
set(MEX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC")
if (APPLE)
  # Trilinos doesn't propagate its sysroot
  message(STATUS "CMAKE_OSX_SYSROOT: ${CMAKE_OSX_SYSROOT}")
  set(MEX_FLAGS "--sysroot=${CMAKE_OSX_SYSROOT} ${MEX_FLAGS}")
endif(APPLE)

list (APPEND CMAKE_INSTALL_RPATH ${Trilinos_LIBRARY_DIRS})
list (APPEND CMAKE_INSTALL_RPATH ${Trilinos_TPL_LIBRARY_DIRS})

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
  # check if library has a "." in it
  string(REGEX MATCH "\.*" lib_ext_result ${alib})
  if (lib_ext_result)

    if (APPLE)
      string(FIND ${alib} ".tbd" alib_result)
      # make sure .tbd not found in library name
    else()
      string(FIND ${alib} ".so" alib_result)
      # make sure .so not found in library name
    endif(APPLE)

    if (alib_result LESS 0)
      # extension found but not .tbd or .so
      message(FATAL_ERROR "${alib} encountered as TPL with unrecognized file extension.")
    else() 
      # remove extension and get library name
      get_filename_component(libname ${alib} NAME_WLE)
      string(SUBSTRING ${libname} 3 -1 reduced_libname)
      set(l_Trilinos_TPL "${l_Trilinos_TPL} -l${reduced_libname}")
    endif()
  else()
    # no extension on library name
    set(l_Trilinos_TPL "${l_Trilinos_TPL} -l${alib}")
  endif()
endforeach()

if(APPLE)
  if(Trilinos_BUILD_SHARED_LIBS)
    set(MEX_LD_FLAGS "LDFLAGS=")
    set(MEX_LD_FLAGS "${MEX_LD_FLAGS} -rpath ${Trilinos_INSTALL_DIR}/lib")
    set(MEX_LD_FLAGS "${MEX_LD_FLAGS} -rpath ${MATLAB_BIN_DIR}/maci64")
    set(MEX_LD_FLAGS "${MEX_LD_FLAGS} -shared")
  endif()
else()
  set(MEX_LD_FLAGS "LDFLAGS=")
  foreach(rp ${CMAKE_INSTALL_RPATH})
    set(MEX_LD_FLAGS "${MEX_LD_FLAGS} -Wl,-R,${rp}")
  endforeach()
  set(MEX_LD_FLAGS "${MEX_LD_FLAGS} -Wl,-R,${MATLAB_BIN_DIR}/glnxa64")
  set(MEX_LD_FLAGS "${MEX_LD_FLAGS} -Wl,--enable-new-dtags")
endif()

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
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/objs
    )

function(make_mex fname)

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o
      COMMAND "${MEX}" "LDFLAGS=-shared" "CXX=${CMAKE_CXX_COMPILER}" "CXXFLAGS=${MEX_FLAGS}" "-largeArrayDims" -c ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp"
      BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/objs
      )
  
  #add_custom_command(OUTPUT ${fname}.${MEX_EXTENSION}
  #    COMMAND "${MEX}" "LDFLAGS=-shared" "CXX=${CMAKE_CXX_COMPILER}" "CXXFLAGS=${MEX_FLAGS}" "-largeArrayDims" ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} ${l_Trilinos} ${l_Trilinos_TPL} "${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o" "${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o"
  #    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o ${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o
  #    BYPRODUCTS ${fname}.${MEX_EXTENSION}
  #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  #    )

  add_custom_command(OUTPUT ${fname}.${MEX_EXTENSION}
      COMMAND "${MEX}" "${MEX_LD_FLAGS}" "CXX=${CMAKE_CXX_COMPILER}" "CXXFLAGS=${MEX_FLAGS}" "-largeArrayDims" ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} ${l_Trilinos} ${l_Trilinos_TPL} "${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/m2i_helpers.cpp"
      #DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/objs/${fname}.o ${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o
      BYPRODUCTS ${fname}.${MEX_EXTENSION}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
  
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${fname}.${MEX_EXTENSION} DESTINATION ${CMAKE_INSTALL_PREFIX})

endfunction(make_mex)

set(MOBJS "intrepid_setJacobian;intrepid_setJacobianInv;intrepid_setJacobianDet;intrepid_getNumCubaturePoints;intrepid_getCubature;intrepid_getBasisValues;intrepid_computeCellMeasure;intrepid_HGRADtransformGRAD;intrepid_HGRADtransformVALUE;intrepid_multiplyMeasure;intrepid_integrate;intrepid_mapToPhysicalFrame;intrepid_mapToReferenceFrame;intrepid_scalarMultiplyDataField;intrepid_scalarMultiplyDataData;intrepid_dataIntegral;intrepid_evaluate;intrepid_multiplyVectorDataField;intrepid_dotMultiplyDataField;intrepid_vectorMultiplyDataField;intrepid_dotMultiplyDataData;intrepid_computeEdgeMeasure;intrepid_computeFaceMeasure;intrepid_getCellDeviatoricMat;intrepid_getCellVolumetricMat;intrepid_getPhysicalSideNormals;intrepid_mapToReferenceSubcell;intrepid_tensorMultiplyDataData;intrepid_tensorMultiplyDataField;defGrad;evaluateVectorField;evaluateVectorGradField;externalForceElastostaticsTL;geometricBmat;intrepid_det;intrepid_inverse;intrepid_transpose;materialBmat;strainEulerianLagrangian;stress;tangentStiffnessElastostaticsTL;tensorMultiplyDataField;myMatVecProductDataField;myMatVecProductDataData")

# dummy targets to build up dependencies
#add_custom_target(m2i_helpers
#    ALL
#    COMMAND # left intentionally blank
#    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/objs/m2i_helpers.o
#    )
foreach(fname ${MOBJS})
    make_mex(${fname})
    # builds dependency chain
    add_custom_target(${fname}
        ALL
        COMMAND # left intentionally blank
        #DEPENDS m2i_helpers ${fname}.${MEX_EXTENSION}
        DEPENDS ${fname}.${MEX_EXTENSION}
        )
endforeach()

if(APPLE)
  message("\"*.dylib was built for newer macOS version (*.*) than being linked (10.9)\" type warnings can be safely ignored. MEX compiler specifies a -mmacosx-min-version=10.9, which can not be overridden.")
endif(APPLE)

