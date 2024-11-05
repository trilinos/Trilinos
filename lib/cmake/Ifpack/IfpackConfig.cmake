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
# CMake variable for use by Trilinos/Ifpack clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Ifpack requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Ifpack build
## ---------------------------------------------------------------------------

set(Ifpack_CXX_COMPILER "/usr/bin/mpicxx")

set(Ifpack_C_COMPILER "/usr/bin/mpicc")

set(Ifpack_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(Ifpack_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Ifpack build
## ---------------------------------------------------------------------------

## Give the build type
set(Ifpack_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Ifpack_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Ifpack_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(Ifpack_Fortran_FLAGS [[ -O3]])
# Deprecated
set(Ifpack_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(Ifpack_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Ifpack_SHARED_LIB_RPATH_COMMAND "")
set(Ifpack_BUILD_SHARED_LIBS "FALSE")

set(Ifpack_LINKER /usr/bin/ld)
set(Ifpack_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Ifpack_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(Ifpack_LIBRARIES Ifpack::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Ifpack_MPI_LIBRARIES "")
set(Ifpack_MPI_LIBRARY_DIRS "")
set(Ifpack_MPI_INCLUDE_DIRS "")
set(Ifpack_MPI_EXEC "/usr/bin/mpiexec")
set(Ifpack_MPI_EXEC_MAX_NUMPROCS "4")
set(Ifpack_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Ifpack_ENABLE_Teuchos ON)
set(Ifpack_ENABLE_Epetra ON)
set(Ifpack_ENABLE_Amesos ON)
set(Ifpack_ENABLE_EpetraExt OFF)
set(Ifpack_ENABLE_AztecOO ON)
set(Ifpack_ENABLE_HYPRE OFF)
set(Ifpack_ENABLE_HIPS OFF)
set(Ifpack_ENABLE_SuperLU OFF)
set(Ifpack_ENABLE_SPARSKIT OFF)
set(Ifpack_ENABLE_Boost OFF)

# Exported cache variables
set(Ifpack_ENABLE_METIS "OFF")
set(HAVE_IFPACK_METIS "OFF")
set(Ifpack_ENABLE_PARALLEL_SUBDOMAIN_SOLVERS "OFF")
set(HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS "OFF")
set(Ifpack_ENABLE_DYNAMIC_FACTORY "OFF")
set(HAVE_IFPACK_DYNAMIC_FACTORY "OFF")
set(Ifpack_ENABLE_Euclid "OFF")
set(HAVE_EUCLID "OFF")
set(Ifpack_ENABLE_TEUCHOS_TIME_MONITOR "ON")
set(IFPACK_TEUCHOS_TIME_MONITOR "ON")
set(Ifpack_ENABLE_FLOPCOUNTERS "OFF")
set(IFPACK_FLOPCOUNTERS "OFF")
set(Ifpack_ENABLE_SUPPORTGRAPH "OFF")
set(HAVE_IFPACK_SUPPORTGRAPH "OFF")
set(Ifpack_ENABLE_Experimental "NO")
set(HAVE_IFPACK_EXPERIMENTAL "OFF")

# Include configuration of dependent packages
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET Amesos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Amesos/AmesosConfig.cmake")
endif()
if (NOT TARGET AztecOO::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../AztecOO/AztecOOConfig.cmake")
endif()

# Import Ifpack targets
include("${CMAKE_CURRENT_LIST_DIR}/IfpackTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Ifpack_IS_TRIBITS_COMPLIANT TRUE)
set(Ifpack_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Ifpack_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Ifpack_EXPORTED_PACKAGE_LIBS_NAMES "ifpack")

foreach(libname IN LISTS Ifpack_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Ifpack::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Ifpack::${libname}', or better yet,"
      " 'Ifpack::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Ifpack'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Ifpack_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
