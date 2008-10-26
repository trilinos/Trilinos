# Including this file will result in the BLAS name mangling to be determined.
#
# NOTE: You should determine the overall Fortran name mangling first before
# you include this file.

IF (F77_FUNC)
  SET(F77_BLAS_MANGLE_DEFAULT ${F77_FUNC})
ELSE()
  SET(F77_BLAS_MANGLE_DEFAULT "UNDEFINED")
ENDIF()

# Set options so that users can change these!

SET(F77_BLAS_MANGLE ${F77_BLAS_MANGLE_DEFAULT} CACHE STRING
  "Name mangling to call functions in the provide BLAS library")

MARK_AS_ADVANCED(F77_BLAS_MANGLE)
