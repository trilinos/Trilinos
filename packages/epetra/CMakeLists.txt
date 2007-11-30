# Project name
PROJECT(Epetra CXX C)

# Setup package version
SET(EPETRA_VERSION "3.7d")

# Flag for turning on MPI support
OPTION(EPETRA_ENABLE_MPI 
       "Build Epetra with MPI for parallel codes."
       OFF
)

# Flag to enable tests
OPTION(EPETRA_ENABLE_TESTS
       "Build Epetra tests."
       OFF
)

# Flag to enable examples
OPTION(EPETRA_ENABLE_EXAMPLES
       "Build Epetra examples."
       OFF
)

# Flag to enable Fortran implementation code for performance
IF(WIN32 AND NOT CYGWIN)
	SET(EPETRA_ENABLE_FORTRAN_DEFAULT OFF)
ELSE(WIN32 AND NOT CYGWIN)
	SET(EPETRA_ENABLE_FORTRAN_DEFAULT ON)
ENDIF(WIN32 AND NOT CYGWIN)

OPTION(EPETRA_ENABLE_FORTRAN
	"Build Fortran library for this package."
	${EPETRA_ENABLE_FORTRAN_DEFAULT}
	)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR/cmake})

# Fortran related configuration
IF(EPETRA_ENABLE_FORTRAN)
  ENABLE_LANGUAGE(Fortran)
ENDIF(EPETRA_ENABLE_FORTRAN)

# MPI related configuration
IF(EPETRA_ENABLE_MPI)
  INCLUDE(ConfigureMPI)
ENDIF(EPETRA_ENABLE_MPI)

# Use secure versions of CRT methods in Windows 
# (addresses C4996 build errors in MSVS)
IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES)
ENDIF(WIN32)

ADD_SUBDIRECTORY(src)

IF(EPETRA_ENABLE_TESTS)
	ADD_SUBDIRECTORY(test)
ENDIF(EPETRA_ENABLE_TESTS)

IF(EPETRA_ENABLE_EXAMPLES)
	ADD_SUBDIRECTORY(example)
ENDIF(EPETRA_ENABLE_EXAMPLES)

