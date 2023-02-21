

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

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#
# Core Files
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  BelosMultiVecTraits_Tpetra.hpp
  BelosOperatorTraits_Tpetra.hpp
  BelosSolverFactory_Tpetra.hpp
  BelosTpetraAdapter.hpp
  BelosTpetraOperator.hpp
  Belos_Details_Tpetra_registerLinearSolverFactory.hpp
  Belos_Details_Tpetra_registerSolverFactory.hpp
  )

APPEND_SET(SOURCES
  BelosTpetraAdapter.cpp
  Belos_Details_Tpetra_registerLinearSolverFactory.cpp
  Belos_Details_Tpetra_registerSolverFactory.cpp
  )

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_Cg.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_CgPipeline.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_CgSingleReduce.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_Gmres.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresS.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresSstep.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresPipeline.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresSingleReduce.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_Krylov.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_Krylov_parameters.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_SolverManager.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_SolverManagerBase.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_UpdateNewton.hpp
  )
APPEND_SET(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_BiCGStab.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_BlockCG.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_BlockGmres.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_Cg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_CgPipeline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_CgSingleReduce.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_FixedPoint.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GCRODR.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_Gmres.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresPipeline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresPoly.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresS.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresSingleReduce.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_GmresSstep.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_LSQR.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_Minres.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_PCPG.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_PseudoBlockCG.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_PseudoBlockGmres.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_PseudoBlockTFQMR.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_TFQMR.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/Belos_Tpetra_RCG.cpp
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
