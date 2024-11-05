# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Including this file will result in the BLAS name mangling to be determined.
#
# NOTE: You should determine the overall Fortran name mangling first before
# you include this file.

if (F77_FUNC)
  set(F77_BLAS_MANGLE_DEFAULT ${F77_FUNC})
else()
  set(F77_BLAS_MANGLE_DEFAULT "UNDEFINED")
endif()
if(WIN32 AND NOT CYGWIN)
  set(F77_BLAS_MANGLE_DEFAULT "${F77_FUNC}")
endif()


# Set options so that users can change these!

set(F77_BLAS_MANGLE ${F77_BLAS_MANGLE_DEFAULT} CACHE STRING
  "Name mangling to call functions in the provided BLAS library")

mark_as_advanced(F77_BLAS_MANGLE)
