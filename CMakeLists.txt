PROJECT(Trilinos)
SET(TRILINOS_VERSION "9.0d")

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Setup user options ...
OPTION(BUILD_SHARED_LIBS "Build shared libraries." OFF)
OPTION(TRILINOS_BUILD_EPETRA "Build the epetra package." ON)
OPTION(TRILINOS_BUILD_TEUCHOS "Build the teuchos package." ON)

# Hide the advanced stuff ...
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES)
MARK_AS_ADVANCED(CMAKE_OSX_SYSROOT)

# Setup global requirements ...
FIND_LIBRARY(TRILINOS_BLAS_LIBRARY NAMES blas blas_win32 DOC "Path to the BLAS implementation")
FIND_LIBRARY(TRILINOS_LAPACK_LIBRARY NAMES lapack lapack_win32 DOC "Path to the LAPACK implementation")

# Setup Fortran name mangling
IF(CYGWIN)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) name ## __")
ENDIF(CYGWIN)

IF(WIN32)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) NAME")
ENDIF(WIN32)

IF(UNIX AND NOT APPLE)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) name ## _")
ENDIF(UNIX AND NOT APPLE)

IF(APPLE)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) name ## __")
ENDIF(APPLE)

# Setup testing ...
INCLUDE(CTest)
MARK_AS_ADVANCED(DART_TESTING_TIMEOUT)

# Setup individual packages ...
ADD_SUBDIRECTORY(packages)

# Setup packaging and distribution ...
SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${TRILINOS_VERSION}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${TRILINOS_VERSION}")
SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${TRILINOS_VERSION}")
SET(CPACK_PACKAGE_NAME "trilinos")
SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
SET(CPACK_PACKAGE_VERSION "${TRILINOS_VERSION}")
SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${TRILINOS_VERSION}")

IF(WIN32)
	SET(CPACK_GENERATOR "NSIS")
ENDIF(WIN32)

INCLUDE(CPack)

