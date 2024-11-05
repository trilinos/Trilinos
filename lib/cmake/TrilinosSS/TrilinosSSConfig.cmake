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
# CMake variable for use by Trilinos/TrilinosSS clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "TrilinosSS requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/TrilinosSS build
## ---------------------------------------------------------------------------

set(TrilinosSS_CXX_COMPILER "/usr/bin/mpicxx")

set(TrilinosSS_C_COMPILER "/usr/bin/mpicc")

set(TrilinosSS_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(TrilinosSS_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/TrilinosSS build
## ---------------------------------------------------------------------------

## Give the build type
set(TrilinosSS_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(TrilinosSS_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(TrilinosSS_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(TrilinosSS_Fortran_FLAGS [[ -O3]])
# Deprecated
set(TrilinosSS_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(TrilinosSS_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(TrilinosSS_SHARED_LIB_RPATH_COMMAND "")
set(TrilinosSS_BUILD_SHARED_LIBS "FALSE")

set(TrilinosSS_LINKER /usr/bin/ld)
set(TrilinosSS_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(TrilinosSS_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(TrilinosSS_LIBRARIES TrilinosSS::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(TrilinosSS_MPI_LIBRARIES "")
set(TrilinosSS_MPI_LIBRARY_DIRS "")
set(TrilinosSS_MPI_INCLUDE_DIRS "")
set(TrilinosSS_MPI_EXEC "/usr/bin/mpiexec")
set(TrilinosSS_MPI_EXEC_MAX_NUMPROCS "4")
set(TrilinosSS_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies

# Exported cache variables
set(TrilinosSS_ENABLE_OpenMP "OFF")
set(TRILINOSSS_HAVE_OMP "OFF")

# Include configuration of dependent packages

# Import TrilinosSS targets
include("${CMAKE_CURRENT_LIST_DIR}/TrilinosSSTargets.cmake")

# Standard TriBITS-compliant external package variables
set(TrilinosSS_IS_TRIBITS_COMPLIANT TRUE)
set(TrilinosSS_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(TrilinosSS_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(TrilinosSS_EXPORTED_PACKAGE_LIBS_NAMES "trilinosss")

foreach(libname IN LISTS TrilinosSS_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE TrilinosSS::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'TrilinosSS::${libname}', or better yet,"
      " 'TrilinosSS::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'TrilinosSS'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'TrilinosSS_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
