# This file gets included in the base-level CMakeLists.txt file to define
# Fortran name mangling.

# Set the default fortran name mangling for each platform

IF(CYGWIN)
  SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
  SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  SET(F77_BLAS_MANGLE_DEFAULT "(name,NAME) name ## _" )
ENDIF()

IF(WIN32)
  SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
  SET(F77_FUNC__DEFAULT "(name,NAME) NAME")
  SET(F77_BLAS__DEFAULTMANGLE "(name,NAME) name ## _" )
ENDIF()

IF(UNIX AND NOT APPLE)
  SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
  #SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  SET(F77_FUNC__DEFAULT "(name,NAME) name ## _" )
  SET(F77_BLAS_MANGLE_DEFAULT "(name,NAME) name ## _" )
ENDIF()

IF(APPLE)
  SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
  SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  SET(F77_BLAS_MANGLE_DEFAULT "(name,NAME) name ## _" )
ENDIF()

# Set options so that users can change these!

SET(F77_FUNC ${F77_FUNC_DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with no underscores in the name")
SET(F77_FUNC_ ${F77_FUNC__DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with at least one underscore in the name")
SET(F77_BLAS_MANGLE ${F77_BLAS_MANGLE_DEFAULT} CACHE STRING
  "Name mangling to call functions in the provide BLAS library")

MARK_AS_ADVANCED(F77_FUNC)
MARK_AS_ADVANCED(F77_FUNC_)
MARK_AS_ADVANCED(F77_BLAS_MANGLE)
