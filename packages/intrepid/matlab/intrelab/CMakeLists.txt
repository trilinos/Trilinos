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


# --> ACTION: define Matlab binary directory
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

# set compilers to Trilinos compilers
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})

# set flags to Trilinos flags
set(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})

# Trilinos doesn't propagate its C++ standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC")
if (APPLE)
  # Trilinos doesn't propagate its sysroot
  message(STATUS "CMAKE_OSX_SYSROOT: ${CMAKE_OSX_SYSROOT}")
  set(CMAKE_CXX_FLAGS "--sysroot=${CMAKE_OSX_SYSROOT} ${CMAKE_CXX_FLAGS}")
endif(APPLE)

# fork of matlab_add_mex from FindMatlab.cmake from CMake@3.21.4
# without symbol exporting and use of *.map files
function(matlab_add_mex_mac)

  if(NOT WIN32)
    # we do not need all this on Windows
    # pthread options
    if(CMAKE_CXX_COMPILER_LOADED)
      check_cxx_compiler_flag(-pthread HAS_MINUS_PTHREAD)
    elseif(CMAKE_C_COMPILER_LOADED)
      check_c_compiler_flag(-pthread HAS_MINUS_PTHREAD)
    endif()
    # we should use try_compile instead, the link flags are discarded from
    # this compiler_flag function.
    #check_cxx_compiler_flag(-Wl,--exclude-libs,ALL HAS_SYMBOL_HIDING_CAPABILITY)

  endif()

  set(options EXECUTABLE MODULE SHARED R2017b R2018a EXCLUDE_FROM_ALL)
  set(oneValueArgs NAME DOCUMENTATION OUTPUT_NAME)
  set(multiValueArgs LINK_TO SRC)

  set(prefix _matlab_addmex_prefix)
  cmake_parse_arguments(${prefix} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT ${prefix}_NAME)
    message(FATAL_ERROR "[MATLAB] The MEX target name cannot be empty")
  endif()

  if(NOT ${prefix}_OUTPUT_NAME)
    set(${prefix}_OUTPUT_NAME ${${prefix}_NAME})
  endif()

  if(NOT Matlab_VERSION_STRING VERSION_LESS "9.1") # For 9.1 (R2016b) and newer, add version source file
    # Add the correct version file depending on which languages are enabled in the project
    if(CMAKE_C_COMPILER_LOADED)
      # If C is enabled, use the .c file as it will work fine also with C++
      set(MEX_VERSION_FILE "${Matlab_ROOT_DIR}/extern/version/c_mexapi_version.c")
    elseif(CMAKE_CXX_COMPILER_LOADED)
      # If C is not enabled, check if CXX is enabled and use the .cpp file
      # to avoid that the .c file is silently ignored
      set(MEX_VERSION_FILE "${Matlab_ROOT_DIR}/extern/version/cpp_mexapi_version.cpp")
    else()
      # If neither C or CXX is enabled, warn because we cannot add the source.
      # TODO: add support for fortran mex files
      message(WARNING "[MATLAB] matlab_add_mex requires that at least C or CXX are enabled languages")
    endif()
  endif()

  # For 9.4 (R2018a) and newer, add API macro.
  # Add it for unknown versions too, just in case.
  if(NOT Matlab_VERSION_STRING VERSION_LESS "9.4"
      OR Matlab_VERSION_STRING STREQUAL "unknown")
    if(${${prefix}_R2018a})
      set(MEX_API_MACRO "MATLAB_DEFAULT_RELEASE=R2018a")
    else()
      set(MEX_API_MACRO "MATLAB_DEFAULT_RELEASE=R2017b")
    endif()
  endif()

  set(_option_EXCLUDE_FROM_ALL)
  if(${prefix}_EXCLUDE_FROM_ALL)
    set(_option_EXCLUDE_FROM_ALL "EXCLUDE_FROM_ALL")
  endif()

  if(${prefix}_EXECUTABLE)
    add_executable(${${prefix}_NAME}
      ${_option_EXCLUDE_FROM_ALL}
      ${${prefix}_SRC}
      ${MEX_VERSION_FILE}
      ${${prefix}_DOCUMENTATION}
      ${${prefix}_UNPARSED_ARGUMENTS})
  else()
    if(${prefix}_MODULE)
      set(type MODULE)
    else()
      set(type SHARED)
    endif()

    add_library(${${prefix}_NAME}
      ${type}
      ${_option_EXCLUDE_FROM_ALL}
      ${${prefix}_SRC}
      ${MEX_VERSION_FILE}
      ${${prefix}_DOCUMENTATION}
      ${${prefix}_UNPARSED_ARGUMENTS})
  endif()

  target_include_directories(${${prefix}_NAME} PRIVATE ${Matlab_INCLUDE_DIRS})

  if(Matlab_HAS_CPP_API)
    if(Matlab_ENGINE_LIBRARY)
      target_link_libraries(${${prefix}_NAME} ${Matlab_ENGINE_LIBRARY})
    endif()
    if(Matlab_DATAARRAY_LIBRARY)
      target_link_libraries(${${prefix}_NAME} ${Matlab_DATAARRAY_LIBRARY})
    endif()
  endif()

  target_link_libraries(${${prefix}_NAME} ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${${prefix}_LINK_TO})
  set_target_properties(${${prefix}_NAME}
      PROPERTIES
        PREFIX ""
        OUTPUT_NAME ${${prefix}_OUTPUT_NAME}
        SUFFIX ".${Matlab_MEX_EXTENSION}")

  target_compile_definitions(${${prefix}_NAME} PRIVATE ${MEX_API_MACRO} MATLAB_MEX_FILE)

endfunction(matlab_add_mex_mac)

set(MOBJS "intrepid_setJacobian;intrepid_setJacobianInv;intrepid_setJacobianDet;intrepid_getNumCubaturePoints;intrepid_getCubature;intrepid_getBasisValues;intrepid_computeCellMeasure;intrepid_HGRADtransformGRAD;intrepid_HGRADtransformVALUE;intrepid_multiplyMeasure;intrepid_integrate;intrepid_mapToPhysicalFrame;intrepid_mapToReferenceFrame;intrepid_scalarMultiplyDataField;intrepid_scalarMultiplyDataData;intrepid_dataIntegral;intrepid_evaluate;intrepid_multiplyVectorDataField;intrepid_dotMultiplyDataField;intrepid_vectorMultiplyDataField;intrepid_dotMultiplyDataData;intrepid_computeEdgeMeasure;intrepid_computeFaceMeasure;intrepid_getCellDeviatoricMat;intrepid_getCellVolumetricMat;intrepid_getPhysicalSideNormals;intrepid_mapToReferenceSubcell;intrepid_tensorMultiplyDataData;intrepid_tensorMultiplyDataField;defGrad;evaluateVectorField;evaluateVectorGradField;externalForceElastostaticsTL;geometricBmat;intrepid_det;intrepid_inverse;intrepid_transpose;materialBmat;strainEulerianLagrangian;stress;tangentStiffnessElastostaticsTL;tensorMultiplyDataField;myMatVecProductDataField;myMatVecProductDataData")

if (APPLE)
  matlab_add_mex_mac(NAME m2i_helpers
      SHARED
      SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/m2i_helpers.cpp
  )
else()
  matlab_add_mex(NAME m2i_helpers
      SHARED
      SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/m2i_helpers.cpp
  )
endif()
target_link_libraries(m2i_helpers Trilinos::all_libs)

function(make_mex fname)
  if (APPLE)
    matlab_add_mex_mac(NAME ${fname}
        SHARED
        SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp
    )
  else()
    matlab_add_mex(NAME ${fname}
        SHARED
        SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/${fname}.cpp
    )
  endif()
  target_link_libraries(${fname} Trilinos::all_libs)
  target_link_libraries(${fname} m2i_helpers)
endfunction(make_mex)

foreach(fname ${MOBJS})
  make_mex(${fname})
endforeach()

install(TARGETS ${MOBJS}
  DESTINATION ${CMAKE_INSTALL_PREFIX} 
)

if(APPLE)
  message("\"*.dylib was built for newer macOS version (*.*) than being linked (10.9)\" type warnings can be safely ignored. MEX compiler specifies a -mmacosx-min-version=10.9, which can not be overridden.")
endif(APPLE)

