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
# CMake variable for use by Trilinos/AztecOO clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "AztecOO requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/AztecOO build
## ---------------------------------------------------------------------------

set(AztecOO_CXX_COMPILER "/usr/bin/mpicxx")

set(AztecOO_C_COMPILER "/usr/bin/mpicc")

set(AztecOO_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(AztecOO_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/AztecOO build
## ---------------------------------------------------------------------------

## Give the build type
set(AztecOO_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(AztecOO_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(AztecOO_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(AztecOO_Fortran_FLAGS [[ -O3]])
# Deprecated
set(AztecOO_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(AztecOO_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(AztecOO_SHARED_LIB_RPATH_COMMAND "")
set(AztecOO_BUILD_SHARED_LIBS "FALSE")

set(AztecOO_LINKER /usr/bin/ld)
set(AztecOO_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(AztecOO_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(AztecOO_LIBRARIES AztecOO::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(AztecOO_MPI_LIBRARIES "")
set(AztecOO_MPI_LIBRARY_DIRS "")
set(AztecOO_MPI_INCLUDE_DIRS "")
set(AztecOO_MPI_EXEC "/usr/bin/mpiexec")
set(AztecOO_MPI_EXEC_MAX_NUMPROCS "4")
set(AztecOO_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(AztecOO_ENABLE_Epetra ON)
set(AztecOO_ENABLE_Teuchos ON)
set(AztecOO_ENABLE_y12m OFF)

# Exported cache variables
set(AztecOO_ENABLE_AZLU "")
set(HAVE_AZLU "OFF")
set(AztecOO_ENABLE_TEUCHOS_TIME_MONITOR "NO")
set(AZ_ENABLE_TIMEMONITOR "OFF")

# Include configuration of dependent packages
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()

# Import AztecOO targets
include("${CMAKE_CURRENT_LIST_DIR}/AztecOOTargets.cmake")

# Standard TriBITS-compliant external package variables
set(AztecOO_IS_TRIBITS_COMPLIANT TRUE)
set(AztecOO_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(AztecOO_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(AztecOO_EXPORTED_PACKAGE_LIBS_NAMES "aztecoo")

foreach(libname IN LISTS AztecOO_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE AztecOO::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'AztecOO::${libname}', or better yet,"
      " 'AztecOO::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'AztecOO'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'AztecOO_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
