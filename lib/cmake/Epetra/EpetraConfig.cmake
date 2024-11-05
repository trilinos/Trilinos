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
# CMake variable for use by Trilinos/Epetra clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Epetra requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Epetra build
## ---------------------------------------------------------------------------

set(Epetra_CXX_COMPILER "/usr/bin/mpicxx")

set(Epetra_C_COMPILER "/usr/bin/mpicc")

set(Epetra_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(Epetra_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Epetra build
## ---------------------------------------------------------------------------

## Give the build type
set(Epetra_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Epetra_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Epetra_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(Epetra_Fortran_FLAGS [[ -O3]])
# Deprecated
set(Epetra_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(Epetra_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Epetra_SHARED_LIB_RPATH_COMMAND "")
set(Epetra_BUILD_SHARED_LIBS "FALSE")

set(Epetra_LINKER /usr/bin/ld)
set(Epetra_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Epetra_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(Epetra_LIBRARIES Epetra::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Epetra_MPI_LIBRARIES "")
set(Epetra_MPI_LIBRARY_DIRS "")
set(Epetra_MPI_INCLUDE_DIRS "")
set(Epetra_MPI_EXEC "/usr/bin/mpiexec")
set(Epetra_MPI_EXEC_MAX_NUMPROCS "4")
set(Epetra_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Epetra_ENABLE_BLAS ON)
set(Epetra_ENABLE_LAPACK ON)
set(Epetra_ENABLE_Teuchos ON)
set(Epetra_ENABLE_CASK OFF)
set(Epetra_ENABLE_MPI ON)
set(Epetra_ENABLE_Oski OFF)

# Exported cache variables
set(Epetra_ENABLE_DEBUG "OFF")
set(HAVE_EPETRA_DEBUG "OFF")
set(Epetra_ENABLE_ABC "OFF")
set(HAVE_EPETRA_ARRAY_BOUNDS_CHECK "OFF")
set(Epetra_ENABLE_FORMAT_IO "OFF")
set(HAVE_FORMAT_IO "OFF")
set(Epetra_ENABLE_WARNING_MESSAGES "OFF")
set(HAVE_WARNING_MESSAGES "OFF")
set(Epetra_ENABLE_FATAL_MESSAGES "OFF")
set(HAVE_FATAL_MESSAGES "OFF")
set(Epetra_ENABLE_THREADS "OFF")
set(HAVE_THREADS "OFF")
set(Epetra_ENABLE_Fortran "ON")
set(HAVE_FORTRAN_SUPPORT "ON")
set(Epetra_DISABLE_READY_SEND_IN_DO_POSTS "OFF")
set(EPETRA_NO_READY_SEND_IN_DO_POSTS "OFF")
set(Trilinos_NO_32BIT_GLOBAL_INDICES "OFF")
set(EPETRA_NO_32BIT_GLOBAL_INDICES "OFF")
set(Trilinos_NO_64BIT_GLOBAL_INDICES "OFF")
set(EPETRA_NO_64BIT_GLOBAL_INDICES "OFF")
set(Epetra_ENABLE_OpenMP "OFF")
set(EPETRA_HAVE_OMP "OFF")

# Include configuration of dependent packages
if (NOT TARGET BLAS::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/BLAS/BLASConfig.cmake")
endif()
if (NOT TARGET LAPACK::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/LAPACK/LAPACKConfig.cmake")
endif()
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET MPI::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/MPI/MPIConfig.cmake")
endif()

# Import Epetra targets
include("${CMAKE_CURRENT_LIST_DIR}/EpetraTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Epetra_IS_TRIBITS_COMPLIANT TRUE)
set(Epetra_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Epetra_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Epetra_EXPORTED_PACKAGE_LIBS_NAMES "epetra")

foreach(libname IN LISTS Epetra_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Epetra::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Epetra::${libname}', or better yet,"
      " 'Epetra::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Epetra'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Epetra_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
