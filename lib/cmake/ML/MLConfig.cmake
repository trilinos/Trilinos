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
# CMake variable for use by Trilinos/ML clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "ML requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/ML build
## ---------------------------------------------------------------------------

set(ML_CXX_COMPILER "/usr/bin/mpicxx")

set(ML_C_COMPILER "/usr/bin/mpicc")

set(ML_Fortran_COMPILER "/usr/bin/mpifort")
# Deprecated!
set(ML_FORTRAN_COMPILER "/usr/bin/mpifort") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/ML build
## ---------------------------------------------------------------------------

## Give the build type
set(ML_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(ML_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(ML_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(ML_Fortran_FLAGS [[ -O3]])
# Deprecated
set(ML_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(ML_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(ML_SHARED_LIB_RPATH_COMMAND "")
set(ML_BUILD_SHARED_LIBS "FALSE")

set(ML_LINKER /usr/bin/ld)
set(ML_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(ML_INSTALL_DIR "/home/as/ThirdParty-v2406/Trilinos")

## List of package libraries
set(ML_LIBRARIES ML::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(ML_MPI_LIBRARIES "")
set(ML_MPI_LIBRARY_DIRS "")
set(ML_MPI_INCLUDE_DIRS "")
set(ML_MPI_EXEC "/usr/bin/mpiexec")
set(ML_MPI_EXEC_MAX_NUMPROCS "4")
set(ML_MPI_EXEC_NUMPROCS_FLAG "-np")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(ML_ENABLE_BLAS ON)
set(ML_ENABLE_LAPACK ON)
set(ML_ENABLE_Teuchos ON)
set(ML_ENABLE_Epetra ON)
set(ML_ENABLE_Zoltan OFF)
set(ML_ENABLE_Galeri OFF)
set(ML_ENABLE_Amesos ON)
set(ML_ENABLE_Ifpack ON)
set(ML_ENABLE_AztecOO ON)
set(ML_ENABLE_EpetraExt OFF)
set(ML_ENABLE_Isorropia OFF)
set(ML_ENABLE_MPI ON)
set(ML_ENABLE_METIS OFF)
set(ML_ENABLE_ParMETIS OFF)
set(ML_ENABLE_PETSC OFF)
set(ML_ENABLE_SuperLU OFF)
set(ML_ENABLE_MATLAB OFF)

# Exported cache variables
set(ML_ENABLE_Aztec "OFF")
set(HAVE_ML_AZTEC "OFF")
set(ML_ENABLE_SuperLU1_0 "OFF")
set(HAVE_ML_SUPERLU1_0 "OFF")
set(ML_ENABLE_SuperLU2_0 "OFF")
set(HAVE_ML_SUPERLU2_0 "OFF")
set(ML_ENABLE_SUPERLUDIST "OFF")
set(HAVE_ML_SUPERLUDIST "OFF")
set(ML_ENABLE_Enrich "OFF")
set(HAVE_ML_ENRICH "OFF")
set(ML_ENABLE_Memory_Check "OFF")
set(HAVE_ML_MEMORY_CHECK "OFF")
set(ML_ENABLE_New_T_PE "OFF")
set(HAVE_ML_NEW_T_PE "OFF")
set(ML_ENABLE_Complex_Maxwell "OFF")
set(HAVE_ML_COMPLEX_MAXWELL "OFF")
set(ML_ENABLE_Timing "OFF")
set(HAVE_ML_TIMING "OFF")
set(ML_ENABLE_Flops "OFF")
set(HAVE_ML_FLOPS "OFF")
set(ML_ENABLE_MLapi "OFF")
set(HAVE_ML_MLAPI "OFF")
set(ML_ENABLE_Cfunc "OFF")
set(HAVE_ML_CFUNC "OFF")
set(ML_ENABLE_TekoSmoothers "OFF")
set(HAVE_ML_TekoSmoothers "OFF")

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
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET Amesos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Amesos/AmesosConfig.cmake")
endif()
if (NOT TARGET Ifpack::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Ifpack/IfpackConfig.cmake")
endif()
if (NOT TARGET AztecOO::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../AztecOO/AztecOOConfig.cmake")
endif()
if (NOT TARGET MPI::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/MPI/MPIConfig.cmake")
endif()

# Import ML targets
include("${CMAKE_CURRENT_LIST_DIR}/MLTargets.cmake")

# Standard TriBITS-compliant external package variables
set(ML_IS_TRIBITS_COMPLIANT TRUE)
set(ML_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(ML_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(ML_EXPORTED_PACKAGE_LIBS_NAMES "ml")

foreach(libname IN LISTS ML_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE ML::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'ML::${libname}', or better yet,"
      " 'ML::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'ML'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'ML_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
