

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
  )

APPEND_SET(SOURCES
  BelosEpetraAdapter.cpp
  BelosEpetraOperator.cpp
  BelosEpetraUtils.cpp
  Belos_Details_Epetra_registerLinearSolverFactory.cpp
  Belos_Details_Epetra_registerSolverFactory.cpp
  Belos_Details_Epetra_registerBiCGStabSolMgr.cpp
  Belos_Details_Epetra_registerBlockCGSolMgr.cpp
  Belos_Details_Epetra_registerBlockGmresSolMgr.cpp
  Belos_Details_Epetra_registerFixedPointSolMgr.cpp
  Belos_Details_Epetra_registerGCRODRSolMgr.cpp
  Belos_Details_Epetra_registerLSQRSolMgr.cpp
  Belos_Details_Epetra_registerPCPGSolMgr.cpp
  Belos_Details_Epetra_registerPseudoBlockCGSolMgr.cpp
  Belos_Details_Epetra_registerPseudoBlockGmresSolMgr.cpp
  Belos_Details_Epetra_registerRCGSolMgr.cpp
  Belos_Details_Epetra_registerMinresSolMgr.cpp
  Belos_Details_Epetra_registerTFQMRSolMgr.cpp
  Belos_Details_Epetra_registerPseudoBlockTFQMRSolMgr.cpp
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
