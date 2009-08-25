# This file gets included in the base-level CMakeLists.txt file to define
# Fortran name mangling.

IF (Trilinos_ENABLE_Fortran)
  INCLUDE(FortranMangling)
  FORTRAN_MANGLING()

  # Verify the selected combination of Fortran and C++ compilers.
  IF("${CMAKE_VERSION}" VERSION_GREATER 2.7.20090824)
    INCLUDE(FortranCInterface)
    FortranCInterface_VERIFY(CXX)
  ENDIF()
ENDIF()

IF (FC_FUNC_DEFAULT)

  SET(F77_FUNC_DEFAULT ${FC_FUNC_DEFAULT})
  SET(F77_FUNC__DEFAULT ${FC_FUNC__DEFAULT})
  # 2008/10/26: rabartl: ToDo: Above, we need to write
  # a different function to find out the right BLAS
  # name mangling automatically.  Given what the above
  # FORTRAN_MANGLING() function does, this should not
  # be too hard.

ELSE()
 
  IF(CYGWIN)
    SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  ELSEIF(WIN32)
    SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    SET(F77_FUNC__DEFAULT "(name,NAME) NAME")
  ELSEIF(UNIX AND NOT APPLE)
    SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    #SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
    SET(F77_FUNC__DEFAULT "(name,NAME) name ## _" )
  ELSEIF(APPLE)
    SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  ELSE()
    MESSAGE(FATAL_ERROR "Error, could not determine fortran name mangling!")
  ENDIF()

ENDIF()

# Set options so that users can change these!

SET(F77_FUNC ${F77_FUNC_DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with no underscores in the name")
SET(F77_FUNC_ ${F77_FUNC__DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with at least one underscore in the name")

MARK_AS_ADVANCED(F77_FUNC)
MARK_AS_ADVANCED(F77_FUNC_)
