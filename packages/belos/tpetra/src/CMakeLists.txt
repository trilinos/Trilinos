
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
  BelosTpetraAdapter.hpp 
  )

APPEND_SET(SOURCES
  BelosTpetraAdapter.cpp
  )

#
# C) Define the targets for package's library(s)
#
# If Kokkos was built with TSQR (Tall Skinny QR) support, link with
# the TSQR libraries.
#
IF (${PACKAGE_NAME}_ENABLE_TSQR})
  PACKAGE_ADD_LIBRARY(
    belostpetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS belos tpetra kokkos kokkosnodetsqr kokkosdisttsqr
    )
ELSE ()
  PACKAGE_ADD_LIBRARY(
    belostpetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS belos tpetra kokkos
    )
ENDIF ()
