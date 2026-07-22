# @HEADER
# *****************************************************************************
#           Trilinos: An Object-Oriented Solver Framework
#
# Copyright 2001-2024 NTESS and the Trilinos contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

if (${CMAKE_VERSION} GREATER "3.13")
     cmake_policy(SET CMP0074 NEW)
endif()

# First, set up the variables for the (backward-compatible) TriBITS way of
# finding ADIOS2.  These are used in case find_package(ADIOS2 ...) is not
# called or does not find HDF5.  Also, these variables need to be non-null
# in order to trigger the right behavior in the function
# tribits_tpl_find_include_dirs_and_libraries().

set(REQUIRED_HEADERS adios2.h)
set(REQUIRED_LIBS_NAMES adios2_c adios2_cxx11)

if (TPL_ENABLE_MPI)
  set(REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} adios2_c_mpi adios2_cxx11_mpi)
endif()

#
# Second, search for ADIOS2 components (if allowed) using the standard
# find_package(ADIOS2 ...).
#
tribits_tpl_allow_pre_find_package(ADIOS2  ADIOS2_ALLOW_PREFIND)
if (ADIOS2_ALLOW_PREFIND)

  message("-- Using find_package(ADIOS2 ...) ...")
  set(CMAKE_MODULE_PATH
    "${CMAKE_MODULE_PATH}"
    "${CMAKE_CURRENT_LIST_DIR}/find_modules"
    "${CMAKE_CURRENT_LIST_DIR}/utils"
     )

  find_package(ADIOS2)
endif()

#
# Third, call tribits_tpl_find_include_dirs_and_libraries()
#
tribits_tpl_find_include_dirs_and_libraries( ADIOS2
  REQUIRED_HEADERS ${REQUIRED_HEADERS}
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )
# NOTE: If find_package(ADIOS2 ...) was called and successfully found ADIOS2,
# then tribits_tpl_find_include_dirs_and_libraries() will use the already-set
# variables TPL_ADIOS2_INCLUDE_DIRS and TPL_ADIOS2_LIBRARIES and then print
# them out (and set some other standard variables as well).  This is the final
# "hook" into the TriBITS TPL system.

# If the `find_package(ADIOS2)` is not run, then this may not be set
# Need to determine how this is set in the library that is being used...

if ("${TPL_ADIOS2_PARALLEL}" STREQUAL "")
   assert_defined(TPL_ADIOS2_INCLUDE_DIRS)
   if ("${TPL_ADIOS2_PARALLEL}" STREQUAL "")
      set(TPL_ADIOS2_PARALLEL False CACHE INTERNAL
          "True if ADIOS2 compiled with parallel enabled")
   endif()
endif()
message(STATUS "TPL_ADIOS2_PARALLEL is ${TPL_ADIOS2_PARALLEL}")
