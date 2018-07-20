

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
  BelosEpetraUtils.h
  Belos_Details_Epetra_registerLinearSolverFactory.hpp
  Belos_Details_Epetra_registerSolverFactory.hpp
  BelosSolverFactory_Epetra.hpp
  )

APPEND_SET(SOURCES
  BelosEpetraAdapter.cpp
  BelosEpetraOperator.cpp
  BelosEpetraUtils.cpp
  Belos_Details_Epetra_registerLinearSolverFactory.cpp
  Belos_Details_Epetra_registerSolverFactory.cpp
  )

#
# C) Define the targets for package's library(s)
#
TRIBITS_ADD_LIBRARY(
  belosepetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS belos
  )
