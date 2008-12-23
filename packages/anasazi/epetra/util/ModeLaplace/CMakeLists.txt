
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

PACKAGE_ADD_LIBRARY(
  ModeLaplace
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS anasazi
  )
