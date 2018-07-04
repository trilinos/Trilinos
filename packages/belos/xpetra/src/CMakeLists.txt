

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
  BelosXpetraAdapter.hpp
  BelosXpetraAdapterOperator.hpp
  BelosSolverFactory_Xpetra.hpp
  BelosXpetraAdapterMultiVector.hpp
  Belos_Details_Xpetra_registerSolverFactory.hpp
  )

APPEND_SET(SOURCES
  Belos_Details_Xpetra_registerSolverFactory.cpp
  )

#
# C) Define the targets for package's library(s)
#
TRIBITS_ADD_LIBRARY(
  belosxpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS belos
  )
