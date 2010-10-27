
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

SET(CURRENT_PACKAGE ANASAZIEPETRA)
CONFIGURE_FILE("${TRILINOS_HOME_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}epetra_DLLExportMacro.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}epetra_DLLExportMacro.h
  )

#
# Core Files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  AnasaziEpetraAdapter.hpp 
  )

APPEND_SET(SOURCES
  AnasaziEpetraAdapter.cpp
  )

#
# C) Define the targets for package's library(s)
#

IF (${PACKAGE_NAME}_ENABLE_TSQR)
  PACKAGE_ADD_LIBRARY(
    anasaziepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS anasazi kokkos kokkosnodetsqr kokkosdisttsqr
    DEFINES -DANASAZIEPETRA_LIB_EXPORTS_MODE
  )
ELSE ()
  PACKAGE_ADD_LIBRARY(
    anasaziepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS anasazi
    DEFINES -DANASAZIEPETRA_LIB_EXPORTS_MODE
    )
ENDIF ()
