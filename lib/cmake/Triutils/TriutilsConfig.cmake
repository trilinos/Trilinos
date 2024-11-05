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
# CMake variable for use by Trilinos/Triutils clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Triutils requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Triutils build
## ---------------------------------------------------------------------------

set(Triutils_CXX_COMPILER "/usr/bin/mpicxx")

set(Triutils_C_COMPILER "/usr/bin/mpicc")

set(Triutils_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(Triutils_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Triutils build
## ---------------------------------------------------------------------------

## Give the build type
set(Triutils_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Triutils_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Triutils_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(Triutils_Fortran_FLAGS [[ -O3]])
# Deprecated
set(Triutils_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(Triutils_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Triutils_SHARED_LIB_RPATH_COMMAND "")
set(Triutils_BUILD_SHARED_LIBS "FALSE")

set(Triutils_LINKER /usr/bin/ld)
set(Triutils_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Triutils_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(Triutils_LIBRARIES Triutils::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Triutils_MPI_LIBRARIES "")
set(Triutils_MPI_LIBRARY_DIRS "")
set(Triutils_MPI_INCLUDE_DIRS "")
set(Triutils_MPI_EXEC "/usr/bin/mpiexec")
set(Triutils_MPI_EXEC_MAX_NUMPROCS "4")
set(Triutils_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Triutils_ENABLE_TeuchosCore ON)
set(Triutils_ENABLE_Epetra ON)

# Exported cache variables

# Include configuration of dependent packages
if (NOT TARGET TeuchosCore::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosCore/TeuchosCoreConfig.cmake")
endif()
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()

# Import Triutils targets
include("${CMAKE_CURRENT_LIST_DIR}/TriutilsTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Triutils_IS_TRIBITS_COMPLIANT TRUE)
set(Triutils_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Triutils_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Triutils_EXPORTED_PACKAGE_LIBS_NAMES "triutils")

foreach(libname IN LISTS Triutils_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Triutils::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Triutils::${libname}', or better yet,"
      " 'Triutils::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Triutils'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Triutils_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
