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
# CMake variable for use by Trilinos/TeuchosParameterList clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "TeuchosParameterList requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/TeuchosParameterList build
## ---------------------------------------------------------------------------

set(TeuchosParameterList_CXX_COMPILER "/usr/bin/mpicxx")

set(TeuchosParameterList_C_COMPILER "/usr/bin/mpicc")

set(TeuchosParameterList_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(TeuchosParameterList_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/TeuchosParameterList build
## ---------------------------------------------------------------------------

## Give the build type
set(TeuchosParameterList_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(TeuchosParameterList_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(TeuchosParameterList_C_FLAGS [[  -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(TeuchosParameterList_Fortran_FLAGS [[ -O3]])
# Deprecated
set(TeuchosParameterList_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(TeuchosParameterList_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(TeuchosParameterList_SHARED_LIB_RPATH_COMMAND "")
set(TeuchosParameterList_BUILD_SHARED_LIBS "FALSE")

set(TeuchosParameterList_LINKER /usr/bin/ld)
set(TeuchosParameterList_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(TeuchosParameterList_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(TeuchosParameterList_LIBRARIES TeuchosParameterList::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(TeuchosParameterList_MPI_LIBRARIES "")
set(TeuchosParameterList_MPI_LIBRARY_DIRS "")
set(TeuchosParameterList_MPI_INCLUDE_DIRS "")
set(TeuchosParameterList_MPI_EXEC "/usr/bin/mpiexec")
set(TeuchosParameterList_MPI_EXEC_MAX_NUMPROCS "4")
set(TeuchosParameterList_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(TeuchosParameterList_ENABLE_TeuchosCore ON)
set(TeuchosParameterList_ENABLE_TeuchosParser ON)
set(TeuchosParameterList_ENABLE_yamlcpp OFF)

# Exported cache variables
set(Teuchos_ENABLE_DEBUG "OFF")
set(HAVE_TEUCHOS_DEBUG "OFF")
set(Teuchos_ENABLE_EXPLICIT_INSTANTIATION "ON")
set(HAVE_TEUCHOS_EXPLICIT_INSTANTIATION "ON")
set(Teuchos_ENABLE_LONG_DOUBLE "FALSE")
set(HAVE_TEUCHOS_LONG_DOUBLE "OFF")
set(Teuchos_ENABLE_FLOAT "FALSE")
set(HAVE_TEUCHOS_FLOAT "OFF")
set(Teuchos_ENABLE_COMPLEX "FALSE")
set(HAVE_TEUCHOS_COMPLEX "OFF")
set(Teuchos_INST_FLOAT "FALSE")
set(HAVE_TEUCHOS_INST_FLOAT "OFF")
set(Teuchos_INST_COMPLEX_FLOAT "OFF")
set(HAVE_TEUCHOS_INST_COMPLEX_FLOAT "OFF")
set(Teuchos_INST_COMPLEX_DOUBLE "FALSE")
set(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE "OFF")
set(Teuchos_ENABLE_THREAD_SAFE "OFF")
set(HAVE_TEUCHOS_THREAD_SAFE "OFF")
set(Teuchos_ENABLE_ABC "OFF")
set(HAVE_TEUCHOS_ARRAY_BOUNDSCHECK "OFF")
set(Teuchos_MODIFY_DEFAULTS_DURING_VALIDATION "ON")
set(HAVE_TEUCHOS_MODIFY_DEFAULTS_DURING_VALIDATION "ON")
set(Teuchos_ENABLE_DEBUG_RCP_NODE_TRACING "OFF")
set(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING "OFF")
set(Teuchos_ENABLE_EXTENDED "ON")
set(HAVE_TEUCHOS_EXTENDED "ON")
set(Teuchos_ENABLE_C_EXCEPTIONS "OFF")
set(HAVE_TEUCHOS_C_EXCEPTIONS "OFF")
set(Teuchos_ENABLE_GCC_DEMANGLE "ON")
set(HAVE_TEUCHOS_DEMANGLE "ON")
set(Teuchos_ENABLE_STACKED_TIMER_IN_TIME_MONITOR "ON")
set(HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER "ON")
set(Teuchos_GLOBALLY_REDUCE_UNITTEST_RESULTS "OFF")
set(HAVE_TEUCHOS_GLOBALLY_REDUCE_UNITTEST_RESULTS "OFF")
set(Teuchos_KOKKOS_PROFILING "ON")
set(HAVE_TEUCHOS_KOKKOS_PROFILING "ON")
set(Teuchos_TIMER_KOKKOS_FENCE "OFF")
set(HAVE_TEUCHOS_TIMER_KOKKOS_FENCE "OFF")

# Include configuration of dependent packages
if (NOT TARGET TeuchosCore::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosCore/TeuchosCoreConfig.cmake")
endif()
if (NOT TARGET TeuchosParser::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosParser/TeuchosParserConfig.cmake")
endif()

# Import TeuchosParameterList targets
include("${CMAKE_CURRENT_LIST_DIR}/TeuchosParameterListTargets.cmake")

# Standard TriBITS-compliant external package variables
set(TeuchosParameterList_IS_TRIBITS_COMPLIANT TRUE)
set(TeuchosParameterList_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(TeuchosParameterList_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(TeuchosParameterList_EXPORTED_PACKAGE_LIBS_NAMES "teuchosparameterlist")

foreach(libname IN LISTS TeuchosParameterList_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE TeuchosParameterList::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'TeuchosParameterList::${libname}', or better yet,"
      " 'TeuchosParameterList::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'TeuchosParameterList'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'TeuchosParameterList_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
