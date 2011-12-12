
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

SET(CURRENT_PACKAGE ANASAZIEPETRA_MODELAPLACE)
CONFIGURE_FILE("${Trilinos_SOURCE_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}epetra_ModeLaplace_DLLExportMacro.h)

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
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}epetra_ModeLaplace_DLLExportMacro.h
  )

#
# Core Files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  BlockPCGSolver.h
  CheckingTools.h
  ModeLaplace1DQ1.h
  ModeLaplace1DQ2.h
  ModeLaplace2DQ1.h
  ModeLaplace2DQ2.h
  ModeLaplace3DQ1.h
  ModeLaplace3DQ2.h
  ModalProblem.h
  ModeLaplace.h
  SortingTools.h
  )

APPEND_SET(SOURCES
  BlockPCGSolver.cpp
  CheckingTools.cpp
  ModeLaplace1DQ1.cpp
  ModeLaplace1DQ2.cpp
  ModeLaplace2DQ1.cpp
  ModeLaplace2DQ2.cpp
  ModeLaplace3DQ1.cpp
  ModeLaplace3DQ2.cpp
  SortingTools.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  ModeLaplace
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS anasazi
  DEFINES -DANASAZIEPETRA_MODELAPLACE_LIB_EXPORTS_MODE
  )
