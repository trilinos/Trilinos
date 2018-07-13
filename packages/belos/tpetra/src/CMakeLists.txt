

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
  BelosMultiVecTraits_Tpetra.hpp
  BelosOperatorTraits_Tpetra.hpp
  BelosSolverFactory_Tpetra.hpp
  BelosTpetraAdapter.hpp
  Belos_Details_Tpetra_registerLinearSolverFactory.hpp
  Belos_Details_Tpetra_registerSolverFactory.hpp
  )

APPEND_SET(SOURCES
  BelosTpetraAdapter.cpp
  Belos_Details_Tpetra_registerLinearSolverFactory.cpp
  Belos_Details_Tpetra_registerSolverFactory.cpp
  )

#
# C) Define the targets for package's library(s)
#
TRIBITS_ADD_LIBRARY(
  belostpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS belos
  )
