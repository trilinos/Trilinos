
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#
# Core Files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  BelosEpetraAdapter.hpp 
  BelosEpetraOperator.h 
  )

APPEND_SET(SOURCES
  BelosEpetraAdapter.cpp
  BelosEpetraOperator.cpp
  )

#
# C) Define the targets for package's library(s)
#

IF (Belos_ENABLE_TSQR)
  PACKAGE_ADD_LIBRARY(
    belosepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS belos epetra kokkos kokkosnodetsqr kokkosdisttsqr
  )
ELSE ()
  PACKAGE_ADD_LIBRARY(
    belosepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS belos epetra
    )
ENDIF ()
