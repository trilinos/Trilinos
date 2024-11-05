# @HEADER
# *****************************************************************************
#           Trilinos: An Object-Oriented Solver Framework
#
# Copyright 2001-2024 NTESS and the Trilinos contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


set(REQUIRED_LIBS_NAMES "blas blas_win32")

#
# Second, search for BLAS components (if allowed) using the standard
# find_package(BLAS ...).
#
tribits_tpl_allow_pre_find_package(BLAS  BLAS_ALLOW_PREFIND)
if (BLAS_ALLOW_PREFIND)

  message("-- Using find_package(BLAS ...) ...")

  find_package(BLAS)

  if (BLAS_FOUND)
    # Tell TriBITS that we found BLAS and there no need to look any further!
    set(TPL_BLAS_INCLUDE_DIRS "" CACHE PATH
      "BLAS include dirs")
    set(TPL_BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE FILEPATH
      "BLAS libraries")
    set(TPL_BLAS_LIBRARY_DIRS "" CACHE PATH
      "BLAS library dirs")
  endif()

endif()

if (MSVC AND NOT
    (BLAS_LIBRARY_DIRS  OR
     (NOT "${BLAS_LIBRARY_NAMES}" STREQUAL "blas blas_win32" AND
      NOT "${BLAS_LIBRARY_NAMES}" STREQUAL "") OR
     BLAS_INCLUDE_DIRS  OR
     BLAS_INCLUDE_NAMES OR
     (NOT "${TPL_BLAS_LIBRARIES}" STREQUAL "blas" AND
      NOT "${TPL_BLAS_LIBRARIES}" STREQUAL "") OR
     TPL_BLAS_INCLUDE_DIRS)
   )
  # Find the CLAPACK built by CMake on the machine for MSVC if found it will
  # set the BLAS and LAPACK libraries.  NOTE: This the FindCLAPACK module must
  # be called every configure or this does not work!
  # If the user has specified alternate name or location of their blas that
  # will be used instead.
  find_package(CLAPACK 3.2.1 NO_MODULE)
  if (CLAPACK_FOUND)
    advanced_set(TPL_BLAS_LIBRARIES blas
      CACHE FILEPATH "Set from MSVC CLAPACK specialization")
  endif()
endif()
#
# Third, call tribits_tpl_find_include_dirs_and_libraries()
#
tribits_tpl_find_include_dirs_and_libraries( BLAS
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )
# NOTE: If find_package(BLAS ...) was called and successfully found BLAS, then
# tribits_tpl_find_include_dirs_and_libraries() will use the already-set
# variables TPL_BLAS_INCLUDE_DIRS and TPL_BLAS_LIBRARIES and then print them
