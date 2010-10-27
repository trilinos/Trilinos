
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

IF (${PACKAGE_NAME}_ENABLE_TSQR)
  # (mfh 27 Oct 2010) TSQR itself lives in Kokkos and does not depend
  # on Tpetra.  However, the Epetra wrappers for TSQR sit in the
  # Tpetra package.  (They do not need to link with the Tpetra
  # library, but we insert a dependency here just to make sure.)  This
  # is a hack to avoid an optional package dependency of Epetra on
  # Teuchos and Kokkos.
  PACKAGE_ADD_LIBRARY(
    belosepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS belos epetra tpetra kokkos kokkosnodetsqr kokkosdisttsqr
  )
ELSE ()
  PACKAGE_ADD_LIBRARY(
    belosepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS belos epetra
    )
ENDIF ()
