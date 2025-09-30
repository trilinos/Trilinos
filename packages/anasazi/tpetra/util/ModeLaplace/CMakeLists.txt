

#
# A) Package-specific configuration options
#

SET(CURRENT_PACKAGE ANASAZITPETRA_MODELAPLACE)
CONFIGURE_FILE("${Trilinos_SOURCE_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}tpetra_ModeLaplace_DLLExportMacro.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}tpetra_ModeLaplace_DLLExportMacro.h
  )

#
# Core Files
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
   TpetraCheckingTools.hpp
   ModeLaplace1DQ1.hpp
   TpetraSortingTools.hpp
  )

APPEND_SET(SOURCES
   TpetraCheckingTools.cpp
   ModeLaplace1DQ1.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  TpetraModeLaplace
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DANASAZITPETRA_MODELAPLACE_LIB_EXPORTS_MODE
  )
