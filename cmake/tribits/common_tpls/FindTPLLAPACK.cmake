# @HEADER
# *****************************************************************************
#           Trilinos: An Object-Oriented Solver Framework
#
# Copyright 2001-2024 NTESS and the Trilinos contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

#
# First, set up the variables for the (backward-compatible) TriBITS way of
# finding LAPACK.  These are used in case find_package(LAPACK ...) is not called
# or does not find LAPACK.  Also, these variables need to be non-null in order
# to trigger the right behavior in the function
# tribits_tpl_find_include_dirs_and_libraries().
#
set(REQUIRED_LIBS_NAMES "lapack lapack_win32")

#
# Second, search for LAPACK components (if allowed) using the standard
# find_package(LAPACK ...).
#
tribits_tpl_allow_pre_find_package(LAPACK  LAPACK_ALLOW_PREFIND)
if (LAPACK_ALLOW_PREFIND)

  message("-- Using find_package(LAPACK ...) ...")

  find_package(LAPACK)

  if (LAPACK_FOUND)
    # Tell TriBITS that we found LAPACK and there no need to look any further!
    set(TPL_LAPACK_INCLUDE_DIRS "" CACHE PATH
      "LAPACK include dirs")
    set(TPL_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE FILEPATH
      "LAPACK libraries")
    set(TPL_LAPACK_LIBRARY_DIRS "" CACHE PATH
      "LAPACK library dirs")
  endif()

endif()

if (MSVC AND NOT
    (LAPACK_LIBRARY_DIRS  OR
     (NOT "${LAPACK_LIBRARY_NAMES}" STREQUAL "lapack lapack_win32" AND
      NOT "${LAPACK_LIBRARY_NAMES}" STREQUAL "") OR
     LAPACK_INCLUDE_DIRS  OR
     LAPACK_INCLUDE_NAMES OR
     (NOT "${TPL_LAPACK_LIBRARIES}" STREQUAL "lapack" AND
      NOT "${TPL_LAPACK_LIBRARIES}" STREQUAL "") OR
     TPL_LAPACK_INCLUDE_DIRS)
   )
  if(CLAPACK_FOUND)
    advanced_set(TPL_LAPACK_LIBRARIES lapack
        CACHE FILEPATH "Set from MSVC CLAPACK specialization")
  endif()
endif()

#
# Third, call tribits_tpl_find_include_dirs_and_libraries()
#
tribits_tpl_find_include_dirs_and_libraries( LAPACK
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )
# NOTE: If find_package(LAPACK ...) was called and successfully found LAPACK, then
# tribits_tpl_find_include_dirs_and_libraries() will use the already-set
# variables TPL_LAPACK_INCLUDE_DIRS and TPL_LAPACK_LIBRARIES and then print them
# out (and set some other standard variables as well).  This is the final
