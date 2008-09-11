# This file gets included in the base-level CMakeLists.txt file to define
# Fortran name mangling.

IF(CYGWIN)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) name ## __")
	SET(F77_BLAS_MANGLE "(name,NAME) name ## _")
ENDIF(CYGWIN)

IF(WIN32)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) NAME")
	SET(F77_BLAS_MANGLE "(name,NAME) name ## _")
ENDIF(WIN32)

IF(UNIX AND NOT APPLE)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) name ## _")
	SET(F77_BLAS_MANGLE "(name,NAME) name ## _")
ENDIF(UNIX AND NOT APPLE)

IF(APPLE)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) name ## __")
	SET(F77_BLAS_MANGLE "(name,NAME) name ## _")
ENDIF(APPLE)

