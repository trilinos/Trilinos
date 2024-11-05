# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# This file gets included in the base-level CMakeLists.txt file to define
# Fortran name mangling.

if (${PROJECT_NAME}_ENABLE_CXX AND ${PROJECT_NAME}_ENABLE_Fortran)
  include(FortranMangling)
  fortran_mangling()

  # Verify the selected combination of Fortran and C++ compilers.
  if(NOT ${PROJECT_NAME}_SKIP_FORTRANCINTERFACE_VERIFY_TEST)
    include(FortranCInterface)
    fortrancinterface_verify(CXX)
  endif()
endif()

if (FC_FUNC_DEFAULT)

  set(F77_FUNC_DEFAULT ${FC_FUNC_DEFAULT})
  set(F77_FUNC__DEFAULT ${FC_FUNC__DEFAULT})
  # 2008/10/26: rabartl: ToDo: Above, we need to write
  # a different function to find out the right BLAS
  # name mangling automatically.  Given what the above
  # fortran_mangling() function does, this should not
  # be too hard.

else()

  if(CYGWIN)
    set(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    set(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  elseif(WIN32)
    set(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    set(F77_FUNC__DEFAULT "(name,NAME) NAME")
  elseif(UNIX AND NOT APPLE)
    set(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    #set(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
    set(F77_FUNC__DEFAULT "(name,NAME) name ## _" )
  elseif(APPLE)
    set(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    set(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  else()
    message(FATAL_ERROR "Error, could not determine fortran name mangling!")
  endif()

endif()

# Set options so that users can change these!

set(F77_FUNC ${F77_FUNC_DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with no underscores in the name")
set(F77_FUNC_ ${F77_FUNC__DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with at least one underscore in the name")

mark_as_advanced(F77_FUNC)
mark_as_advanced(F77_FUNC_)
