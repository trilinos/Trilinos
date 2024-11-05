# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

##############################################################################
#
# CMake variable for use by Trilinos/Kokkos clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Kokkos requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Kokkos build
## ---------------------------------------------------------------------------

set(Kokkos_CXX_COMPILER "/usr/bin/mpicxx")

set(Kokkos_C_COMPILER "/usr/bin/mpicc")

set(Kokkos_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(Kokkos_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Kokkos build
## ---------------------------------------------------------------------------

## Give the build type
set(Kokkos_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Kokkos_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Kokkos_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(Kokkos_Fortran_FLAGS [[ -O3]])
# Deprecated
set(Kokkos_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(Kokkos_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Kokkos_SHARED_LIB_RPATH_COMMAND "")
set(Kokkos_BUILD_SHARED_LIBS "FALSE")

set(Kokkos_LINKER /usr/bin/ld)
set(Kokkos_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Kokkos_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(Kokkos_LIBRARIES Kokkos::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Kokkos_MPI_LIBRARIES "")
set(Kokkos_MPI_LIBRARY_DIRS "")
set(Kokkos_MPI_INCLUDE_DIRS "")
set(Kokkos_MPI_EXEC "/usr/bin/mpiexec")
set(Kokkos_MPI_EXEC_MAX_NUMPROCS "4")
set(Kokkos_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Kokkos_ENABLE_Pthread OFF)
set(Kokkos_ENABLE_CUDA OFF)
set(Kokkos_ENABLE_HWLOC OFF)
set(Kokkos_ENABLE_DLlib ON)

# Exported cache variables
set(Kokkos_CXX_STANDARD "")
set(Kokkos_ENABLE_THREADS "OFF")
set(Kokkos_ENABLE_OPENMP "OFF")
set(Kokkos_ENABLE_SERIAL "ON")
set(Kokkos_ENABLE_HPX "OFF")
set(Kokkos_ENABLE_OPENACC "OFF")
set(Kokkos_ENABLE_OPENMPTARGET "OFF")
set(Kokkos_ENABLE_CUDA "")
set(Kokkos_ENABLE_HIP "OFF")
set(Kokkos_ENABLE_SYCL "OFF")
set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE "OFF")
set(Kokkos_ENABLE_CUDA_UVM "OFF")
set(Kokkos_ENABLE_CUDA_LDG_INTRINSIC "OFF")
set(Kokkos_ENABLE_CUDA_LAMBDA "OFF")
set(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC "ON")
set(Kokkos_ENABLE_IMPL_NVHPC_AS_DEVICE_COMPILER "OFF")
set(Kokkos_ENABLE_IMPL_CUDA_UNIFIED_MEMORY "OFF")
set(Kokkos_ENABLE_DEPRECATED_CODE_4 "ON")
set(Kokkos_ENABLE_DEPRECATION_WARNINGS "ON")
set(Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE "OFF")
set(Kokkos_ENABLE_DEBUG "OFF")
set(Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK "OFF")
set(Kokkos_ENABLE_LARGE_MEM_TESTS "OFF")
set(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK "OFF")
set(Kokkos_ENABLE_TUNING "OFF")
set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION "OFF")
set(Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE "OFF")
set(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS "OFF")
set(Kokkos_ENABLE_IMPL_HIP_UNIFIED_MEMORY "OFF")
set(Kokkos_ENABLE_DESUL_ATOMICS_EXTERNAL "OFF")
set(Kokkos_ENABLE_ATOMICS_BYPASS "OFF")
set(Kokkos_ENABLE_IMPL_REF_COUNT_BRANCH_UNLIKELY "ON")
set(Kokkos_ENABLE_IMPL_VIEW_OF_VIEWS_DESTRUCTOR_PRECONDITION_VIOLATION_WORKAROUND "OFF")
set(Kokkos_ENABLE_IMPL_MDSPAN "ON")
set(Kokkos_ENABLE_MDSPAN_EXTERNAL "OFF")
set(Kokkos_ENABLE_IMPL_SKIP_COMPILER_MDSPAN "ON")
set(Kokkos_ENABLE_COMPLEX_ALIGN "OFF")
set(Kokkos_ENABLE_CUDA_CONSTEXPR "OFF")
set(Kokkos_ENABLE_IMPL_HPX_ASYNC_DISPATCH "OFF")
set(Kokkos_ENABLE_UNSUPPORTED_ARCHS "OFF")
set(Kokkos_IMPL_AMDGPU_FLAGS "")
set(Kokkos_IMPL_AMDGPU_LINK "")
set(Kokkos_ENABLE_HWLOC "")
set(Kokkos_HWLOC_DIR "")
set(Kokkos_ENABLE_CUDA "")
set(Kokkos_CUDA_DIR "")
set(Kokkos_ENABLE_ROCM "OFF")
set(Kokkos_ROCM_DIR "")
set(Kokkos_ENABLE_ROCTHRUST "OFF")
set(Kokkos_ROCTHRUST_DIR "")
set(Kokkos_ENABLE_ONEDPL "OFF")
set(Kokkos_ONEDPL_DIR "")
set(Kokkos_ENABLE_LIBDL "ON")
set(Kokkos_LIBDL_DIR "")
set(Kokkos_ENABLE_HPX "OFF")
set(Kokkos_HPX_DIR "")
set(Kokkos_ENABLE_THREADS "OFF")
set(Kokkos_THREADS_DIR "")
set(Kokkos_ENABLE_LIBQUADMATH "OFF")
set(Kokkos_LIBQUADMATH_DIR "")

# Include configuration of dependent packages
if (NOT TARGET DLlib::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/DLlib/DLlibConfig.cmake")
endif()

# Import Kokkos targets
include("${CMAKE_CURRENT_LIST_DIR}/KokkosTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Kokkos_IS_TRIBITS_COMPLIANT TRUE)
set(Kokkos_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Kokkos_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Kokkos_EXPORTED_PACKAGE_LIBS_NAMES "kokkoscore;kokkoscontainers;kokkosalgorithms;kokkossimd")

foreach(libname IN LISTS Kokkos_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Kokkos::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Kokkos::${libname}', or better yet,"
      " 'Kokkos::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Kokkos'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Kokkos_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
SET(Kokkos_DEVICES SERIAL)
SET(Kokkos_OPTIONS IMPL_CUDA_MALLOC_ASYNC;DEPRECATED_CODE_4;DEPRECATION_WARNINGS;IMPL_REF_COUNT_BRANCH_UNLIKELY;IMPL_MDSPAN;IMPL_SKIP_COMPILER_MDSPAN)
SET(Kokkos_TPLS )
SET(Kokkos_ARCH )
SET(Kokkos_CXX_COMPILER "/usr/bin/mpicxx")
SET(Kokkos_CXX_COMPILER_ID "GNU")
SET(Kokkos_CXX_COMPILER_VERSION "13.2.0")
SET(Kokkos_CXX_STANDARD 17)

# Required to be a TriBITS-compliant external package
IF(NOT TARGET Kokkos::all_libs)
  # CMake Error at <prefix>/lib/cmake/Kokkos/KokkosConfigCommon.cmake:10 (ADD_LIBRARY):
  #   ADD_LIBRARY cannot create ALIAS target "Kokkos::all_libs" because target
  #   "Kokkos::kokkos" is imported but not globally visible.
  IF(CMAKE_VERSION VERSION_LESS "3.18")
    SET_TARGET_PROPERTIES(Kokkos::kokkos PROPERTIES IMPORTED_GLOBAL ON)
  ENDIF()
  ADD_LIBRARY(Kokkos::all_libs ALIAS Kokkos::kokkos)
ENDIF()

# Export Kokkos_ENABLE_<BACKEND> for each backend that was enabled.
# NOTE: "Devices" is a little bit of a misnomer here.  These are really
# backends, e.g. Kokkos_ENABLE_OPENMP, Kokkos_ENABLE_CUDA, Kokkos_ENABLE_HIP,
# or Kokkos_ENABLE_SYCL.
FOREACH(DEV ${Kokkos_DEVICES})
  SET(Kokkos_ENABLE_${DEV} ON)
ENDFOREACH()
# Export relevant Kokkos_ENABLE<OPTION> variables, e.g.
# Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE, Kokkos_ENABLE_DEBUG, etc.
FOREACH(OPT ${Kokkos_OPTIONS})
  SET(Kokkos_ENABLE_${OPT} ON)
ENDFOREACH()

IF(Kokkos_ENABLE_CUDA)
  SET(Kokkos_CUDA_ARCHITECTURES )
ENDIF()

IF(Kokkos_ENABLE_HIP)
  SET(Kokkos_HIP_ARCHITECTURES )
ENDIF()

IF(NOT Kokkos_FIND_QUIETLY)
  MESSAGE(STATUS "Enabled Kokkos devices: ${Kokkos_DEVICES}")
ENDIF()

IF (Kokkos_ENABLE_CUDA)
  # If we are building CUDA, we have tricked CMake because we declare a CXX project
  # If the default C++ standard for a given compiler matches the requested
  # standard, then CMake just omits the -std flag in later versions of CMake
  # This breaks CUDA compilation (CUDA compiler can have a different default
  # -std then the underlying host compiler by itself). Setting this variable
  # forces CMake to always add the -std flag even if it thinks it doesn't need it
  SET(CMAKE_CXX_STANDARD_DEFAULT 98 CACHE INTERNAL "" FORCE)
ENDIF()

SET(KOKKOS_USE_CXX_EXTENSIONS OFF)
IF (NOT DEFINED CMAKE_CXX_EXTENSIONS OR CMAKE_CXX_EXTENSIONS)
  IF (NOT KOKKOS_USE_CXX_EXTENSIONS)
    MESSAGE(WARNING "The installed Kokkos configuration does not support CXX extensions. Forcing -DCMAKE_CXX_EXTENSIONS=Off")
    SET(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "" FORCE)
  ENDIF()
ENDIF()

include(FindPackageHandleStandardArgs)

#   This function makes sure that Kokkos was built with the requested backends
#   and target architectures and generates a fatal error if it was not.
#
#   kokkos_check(
#     [DEVICES <devices>...]   # Set of backends (e.g. "OpenMP" and/or "Cuda")
#     [ARCH <archs>...]        # Target architectures (e.g. "Power9" and/or "Volta70")
#     [OPTIONS <options>...]   # Optional settings (e.g. "TUNING")
#     [TPLS <tpls>...]         # Third party libraries
#     [RETURN_VALUE <result>]  # Set a variable that indicates the result of the
#                              # check instead of a fatal error
#   )
function(kokkos_check)
  set(ALLOWED_ARGS DEVICES ARCH OPTIONS TPLS)
  cmake_parse_arguments(KOKKOS_CHECK "" "RETURN_VALUE" "${ALLOWED_ARGS}" ${ARGN})
  foreach(_arg ${KOKKOS_CHECK_UNPARSED_ARGUMENTS})
    message(SEND_ERROR "Argument '${_arg}' passed to kokkos_check() was not recognized")
  endforeach()
  # Get the list of keywords that were actually passed to the function.
  set(REQUESTED_ARGS)
  foreach(arg ${ALLOWED_ARGS})
    if(KOKKOS_CHECK_${arg})
      list(APPEND REQUESTED_ARGS ${arg})
    endif()
  endforeach()
  set(KOKKOS_CHECK_SUCCESS TRUE)
  foreach(arg ${REQUESTED_ARGS})
    # Define variables named after the required arguments that are provided by
    # the Kokkos install.
    foreach(requested ${KOKKOS_CHECK_${arg}})
      foreach(provided ${Kokkos_${arg}})
        STRING(TOUPPER ${requested} REQUESTED_UC)
        STRING(TOUPPER ${provided}  PROVIDED_UC)
        if(PROVIDED_UC STREQUAL REQUESTED_UC)
          string(REPLACE ";" " " ${requested} "${KOKKOS_CHECK_${arg}}")
        endif()
      endforeach()
    endforeach()
    # Somewhat divert the CMake function below from its original purpose and
    # use it to check that there are variables defined for all required
    # arguments. Success or failure messages will be displayed but we are
    # responsible for signaling failure and skip the build system generation.
    if (KOKKOS_CHECK_RETURN_VALUE)
      set(Kokkos_${arg}_FIND_QUIETLY ON)
    endif()
    find_package_handle_standard_args("Kokkos_${arg}" DEFAULT_MSG
            ${KOKKOS_CHECK_${arg}})
    if(NOT Kokkos_${arg}_FOUND)
      set(KOKKOS_CHECK_SUCCESS FALSE)
    endif()
  endforeach()
  if(NOT KOKKOS_CHECK_SUCCESS AND NOT KOKKOS_CHECK_RETURN_VALUE)
    message(FATAL_ERROR "Kokkos does NOT provide all backends and/or architectures requested")
  else()
    set(${KOKKOS_CHECK_RETURN_VALUE} ${KOKKOS_CHECK_SUCCESS} PARENT_SCOPE)
  endif()
endfunction()

# A test to check whether a downstream project set the C++ compiler to NVCC or not
# this is called only when Kokkos was installed with Kokkos_ENABLE_CUDA=ON
FUNCTION(kokkos_compiler_is_nvcc VAR COMPILER)
    # Check if the compiler is nvcc (which really means nvcc_wrapper).
    EXECUTE_PROCESS(COMMAND ${COMPILER} ${ARGN} --version
                    OUTPUT_VARIABLE INTERNAL_COMPILER_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE RET)
    # something went wrong
    IF(RET GREATER 0)
        SET(${VAR} false PARENT_SCOPE)
    ELSE()
        STRING(REPLACE "\n" " - " INTERNAL_COMPILER_VERSION_ONE_LINE ${INTERNAL_COMPILER_VERSION} )
        STRING(FIND ${INTERNAL_COMPILER_VERSION_ONE_LINE} "nvcc" INTERNAL_COMPILER_VERSION_CONTAINS_NVCC)
        STRING(REGEX REPLACE "^ +" "" INTERNAL_HAVE_COMPILER_NVCC "${INTERNAL_HAVE_COMPILER_NVCC}")
        IF(${INTERNAL_COMPILER_VERSION_CONTAINS_NVCC} GREATER -1)
            SET(${VAR} true PARENT_SCOPE)
        ELSE()
            SET(${VAR} false PARENT_SCOPE)
        ENDIF()
    ENDIF()
ENDFUNCTION()

# this function checks whether the current CXX compiler supports building CUDA
FUNCTION(kokkos_cxx_compiler_cuda_test _VAR _COMPILER)

    FILE(WRITE ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cu
"
#include <cuda.h>
#include <cstdlib>

__global__
void kernel(int sz, double* data)
{
    int _beg = blockIdx.x * blockDim.x + threadIdx.x;
    for(int i = _beg; i < sz; ++i)
        data[i] += static_cast<double>(i);
}

int main()
{
    double* data = NULL;
    int blocks = 64;
    int grids = 64;
    int ret = cudaMalloc(&data, blocks * grids * sizeof(double));
    if(ret != cudaSuccess)
        return EXIT_FAILURE;
    kernel<<<grids, blocks>>>(blocks * grids, data);
    cudaDeviceSynchronize();
    return EXIT_SUCCESS;
}
")

    # save the command for debugging
    SET(_COMMANDS "${_COMPILER} ${ARGN} -c ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cu")

    # use execute_process instead of try compile because we want to set custom compiler
    EXECUTE_PROCESS(COMMAND ${_COMPILER} ${ARGN} -c ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cu
        RESULT_VARIABLE     _RET
        WORKING_DIRECTORY   ${PROJECT_BINARY_DIR}/compile_tests
        TIMEOUT             15
        OUTPUT_QUIET
        ERROR_QUIET)

    IF(NOT _RET EQUAL 0)
        # save the command for debugging
        SET(_COMMANDS "${_COMMAND}\n${_COMPILER} --cuda-gpu-arch=sm_35 ${ARGN} -c ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cu")
        # try the compile test again with clang arguments
        EXECUTE_PROCESS(COMMAND ${_COMPILER} --cuda-gpu-arch=sm_35 -c ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cu
            RESULT_VARIABLE     _RET
            WORKING_DIRECTORY   ${PROJECT_BINARY_DIR}/compile_tests
            TIMEOUT             15
            OUTPUT_QUIET
            ERROR_QUIET)
    ENDIF()

    SET(${_VAR}_COMMANDS "${_COMMANDS}" PARENT_SCOPE)
    SET(${_VAR} ${_RET} PARENT_SCOPE)
ENDFUNCTION()

# this function is provided to easily select which files use the same compiler as Kokkos
# when it was installed (or nvcc_wrapper):
#
#       GLOBAL      --> all files
#       TARGET      --> all files in a target
#       SOURCE      --> specific source files
#       DIRECTORY   --> all files in directory
#       PROJECT     --> all files/targets in a project/subproject
#
# Use the COMPILER argument to specify a compiler, if needed. By default, it will
# set the values to ${Kokkos_CXX_COMPILER} unless Kokkos_ENABLE_CUDA=ON and
# Kokkos_CXX_COMPILER_ID is NVIDIA, then it will set it to nvcc_wrapper
#
# Use CHECK_CUDA_COMPILES to run a check when CUDA is enabled
#
FUNCTION(kokkos_compilation)
    CMAKE_PARSE_ARGUMENTS(COMP
        "GLOBAL;PROJECT;CHECK_CUDA_COMPILES"
        "COMPILER"
        "DIRECTORY;TARGET;SOURCE;COMMAND_PREFIX"
        ${ARGN})

    # if built w/o CUDA support, we want to basically make this a no-op
    SET(_Kokkos_ENABLE_CUDA )


    IF(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
      SET(MAYBE_CURRENT_INSTALLATION_ROOT "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../../..")
    ENDIF()

    # search relative first and then absolute
    SET(_HINTS "${MAYBE_CURRENT_INSTALLATION_ROOT}" "/home/as/ThirdParty-v2406/Trilinos")

    # find kokkos_launch_compiler
    FIND_PROGRAM(Kokkos_COMPILE_LAUNCHER
        NAMES           kokkos_launch_compiler
        HINTS           ${_HINTS}
        PATHS           ${_HINTS}
        PATH_SUFFIXES   bin)

    IF(NOT Kokkos_COMPILE_LAUNCHER)
        MESSAGE(FATAL_ERROR "Kokkos could not find 'kokkos_launch_compiler'. Please set '-DKokkos_COMPILE_LAUNCHER=/path/to/launcher'")
    ENDIF()

    # if COMPILER was not specified, assume Kokkos_CXX_COMPILER
    IF(NOT COMP_COMPILER)
        SET(COMP_COMPILER ${Kokkos_CXX_COMPILER})
        IF(_Kokkos_ENABLE_CUDA AND Kokkos_CXX_COMPILER_ID STREQUAL NVIDIA)
            # find nvcc_wrapper
            FIND_PROGRAM(Kokkos_NVCC_WRAPPER
                NAMES           nvcc_wrapper
                HINTS           ${_HINTS}
                PATHS           ${_HINTS}
                PATH_SUFFIXES   bin)
            # fatal if we can't nvcc_wrapper
            IF(NOT Kokkos_NVCC_WRAPPER)
                MESSAGE(FATAL_ERROR "Kokkos could not find nvcc_wrapper. Please set '-DKokkos_NVCC_WRAPPER=/path/to/nvcc_wrapper'")
            ENDIF()
            SET(COMP_COMPILER ${Kokkos_NVCC_WRAPPER})
        ENDIF()
    ENDIF()

    # check that the original compiler still exists!
    IF(NOT EXISTS ${COMP_COMPILER})
        MESSAGE(FATAL_ERROR "Kokkos could not find original compiler: '${COMP_COMPILER}'")
    ENDIF()

    # try to ensure that compiling cuda code works!
    IF(_Kokkos_ENABLE_CUDA AND COMP_CHECK_CUDA_COMPILES)

        # this may fail if kokkos_compiler launcher was used during install
        kokkos_cxx_compiler_cuda_test(_COMPILES_CUDA
            ${Kokkos_COMPILE_LAUNCHER} ${COMP_COMPILER} ${CMAKE_CXX_COMPILER})

        # if above failed, throw an error
        IF(NOT _COMPILES_CUDA)
            MESSAGE(FATAL_ERROR "kokkos_cxx_compiler_cuda_test failed! Test commands:\n${_COMPILES_CUDA_COMMANDS}")
        ENDIF()
    ENDIF()

    IF(COMP_COMMAND_PREFIX)
        SET(_PREFIX "${COMP_COMMAND_PREFIX}")
        STRING(REPLACE ";" " " _PREFIX "${COMP_COMMAND_PREFIX}")
        SET(Kokkos_COMPILER_LAUNCHER "${_PREFIX} ${Kokkos_COMPILE_LAUNCHER}")
    ENDIF()

    IF(COMP_GLOBAL)
        # if global, don't bother setting others
        SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${Kokkos_COMPILE_LAUNCHER} ${COMP_COMPILER} ${CMAKE_CXX_COMPILER}")
        SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_LINK "${Kokkos_COMPILE_LAUNCHER} ${COMP_COMPILER} ${CMAKE_CXX_COMPILER}")
    ELSE()
        FOREACH(_TYPE PROJECT DIRECTORY TARGET SOURCE)
            # make project/subproject scoping easy, e.g. KokkosCompilation(PROJECT) after project(...)
            IF("${_TYPE}" STREQUAL "PROJECT" AND COMP_${_TYPE})
                LIST(APPEND COMP_DIRECTORY ${PROJECT_SOURCE_DIR})
                UNSET(COMP_${_TYPE})
            ENDIF()
            # set the properties if defined
            IF(COMP_${_TYPE})
                # MESSAGE(STATUS "Using ${COMP_COMPILER} :: ${_TYPE} :: ${COMP_${_TYPE}}")
                SET_PROPERTY(${_TYPE} ${COMP_${_TYPE}} PROPERTY RULE_LAUNCH_COMPILE "${Kokkos_COMPILE_LAUNCHER} ${COMP_COMPILER} ${CMAKE_CXX_COMPILER}")
                SET_PROPERTY(${_TYPE} ${COMP_${_TYPE}} PROPERTY RULE_LAUNCH_LINK "${Kokkos_COMPILE_LAUNCHER} ${COMP_COMPILER} ${CMAKE_CXX_COMPILER}")
            ENDIF()
        ENDFOREACH()
    ENDIF()
ENDFUNCTION()
IF (NOT TARGET Kokkos::kokkos)
  # Compute the installation prefix relative to this file.
  get_filename_component(KOKKOS_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
  get_filename_component(KOKKOS_IMPORT_PREFIX "${KOKKOS_IMPORT_PREFIX}" PATH)
  get_filename_component(KOKKOS_IMPORT_PREFIX "${KOKKOS_IMPORT_PREFIX}" PATH)
  get_filename_component(KOKKOS_IMPORT_PREFIX "${KOKKOS_IMPORT_PREFIX}" PATH)
  if(KOKKOS_IMPORT_PREFIX STREQUAL "/")
    set(KOKKOS_IMPORT_PREFIX "")
  endif()
  add_library(Kokkos::kokkos INTERFACE IMPORTED)
  set_target_properties(Kokkos::kokkos PROPERTIES
    INTERFACE_LINK_LIBRARIES "Kokkos::kokkossimd;Kokkos::kokkosalgorithms;Kokkos::kokkoscontainers;Kokkos::kokkoscore;-DKOKKOS_DEPENDENCE"
    INTERFACE_COMPILE_FEATURES "cxx_std_17"
    INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CXX>:>"
    INTERFACE_INCLUDE_DIRECTORIES "${KOKKOS_IMPORT_PREFIX}/include"
  )
ENDIF()
