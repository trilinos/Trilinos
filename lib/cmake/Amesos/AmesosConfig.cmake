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
# CMake variable for use by Trilinos/Amesos clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Amesos requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Amesos build
## ---------------------------------------------------------------------------

set(Amesos_CXX_COMPILER "/usr/bin/mpicxx")

set(Amesos_C_COMPILER "/usr/bin/mpicc")

set(Amesos_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(Amesos_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Amesos build
## ---------------------------------------------------------------------------

## Give the build type
set(Amesos_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Amesos_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Amesos_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(Amesos_Fortran_FLAGS [[ -O3]])
# Deprecated
set(Amesos_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(Amesos_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Amesos_SHARED_LIB_RPATH_COMMAND "")
set(Amesos_BUILD_SHARED_LIBS "FALSE")

set(Amesos_LINKER /usr/bin/ld)
set(Amesos_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Amesos_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(Amesos_LIBRARIES Amesos::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Amesos_MPI_LIBRARIES "")
set(Amesos_MPI_LIBRARY_DIRS "")
set(Amesos_MPI_INCLUDE_DIRS "")
set(Amesos_MPI_EXEC "/usr/bin/mpiexec")
set(Amesos_MPI_EXEC_MAX_NUMPROCS "4")
set(Amesos_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Amesos_ENABLE_Teuchos ON)
set(Amesos_ENABLE_Epetra ON)
set(Amesos_ENABLE_TrilinosSS ON)
set(Amesos_ENABLE_EpetraExt OFF)
set(Amesos_ENABLE_SuperLUDist OFF)
set(Amesos_ENABLE_ParMETIS OFF)
set(Amesos_ENABLE_UMFPACK OFF)
set(Amesos_ENABLE_SuperLU OFF)
set(Amesos_ENABLE_BLACS OFF)
set(Amesos_ENABLE_SCALAPACK OFF)
set(Amesos_ENABLE_MUMPS OFF)
set(Amesos_ENABLE_TAUCS OFF)
set(Amesos_ENABLE_CSS_MKL OFF)
set(Amesos_ENABLE_PARDISO_MKL OFF)
set(Amesos_ENABLE_PARDISO OFF)
set(Amesos_ENABLE_CSparse OFF)

# Exported cache variables
set(Amesos_ENABLE_DSCPACK "OFF")
set(HAVE_AMESOS_DSCPACK "OFF")
set(Amesos_ENABLE_MC64 "OFF")
set(HAVE_AMESOS_MC64 "OFF")
set(Amesos_ENABLE_KLU "ON")
set(HAVE_AMESOS_KLU "ON")
set(Amesos_ENABLE_PARAKLETE "OFF")
set(HAVE_AMESOS_PARAKLETE "OFF")
set(Amesos_ENABLE_LAPACK "ON")
set(HAVE_AMESOS_LAPACK "ON")

# Include configuration of dependent packages
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET TrilinosSS::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TrilinosSS/TrilinosSSConfig.cmake")
endif()

# Import Amesos targets
include("${CMAKE_CURRENT_LIST_DIR}/AmesosTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Amesos_IS_TRIBITS_COMPLIANT TRUE)
set(Amesos_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Amesos_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Amesos_EXPORTED_PACKAGE_LIBS_NAMES "amesos")

foreach(libname IN LISTS Amesos_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Amesos::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Amesos::${libname}', or better yet,"
      " 'Amesos::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Amesos'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Amesos_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
