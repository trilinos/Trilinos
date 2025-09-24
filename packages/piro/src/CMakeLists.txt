
#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#


TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

#
# B.1) utilities
#
APPEND_SET(HEADERS
  Piro_Helpers.hpp 
  Piro_ValidPiroParameters.hpp
  Piro_Version.hpp
  Piro_InvertMassMatrixDecorator.hpp
  Piro_InvertMassMatrixDecorator_Def.hpp
  Piro_Provider.hpp
  Piro_ProviderBase.hpp
  Piro_ProviderHelpers.hpp
  Piro_SolutionObserverBase.hpp
  Piro_SolverFactory.hpp
  Piro_SolverFactory_Def.hpp
  Piro_StratimikosUtils.hpp
  Piro_ConfigDefs.hpp
  )
IF(PIRO_HAVE_EPETRA_STACK)
APPEND_SET(HEADERS
  Piro_Epetra_MatrixFreeDecorator.hpp
  Piro_Epetra_InvertMassMatrixDecorator.hpp
  Piro_Epetra_SolverFactory.hpp
  )
ENDIF()
APPEND_SET(SOURCES
  Piro_ValidPiroParameters.cpp
  Piro_StratimikosUtils.cpp
  Piro_Version.cpp
  )
IF(PIRO_HAVE_EPETRA_STACK)
APPEND_SET(SOURCES
  Piro_Epetra_MatrixFreeDecorator.cpp
  Piro_Epetra_InvertMassMatrixDecorator.cpp
  Piro_Epetra_SolverFactory.cpp
  )
ENDIF()
APPEND_SET(HEADERS
  Piro_MatrixFreeDecorator.hpp
  Piro_MatrixFreeDecorator_Def.hpp
  Piro_MatrixFreeLinearOp.hpp
  Piro_MatrixFreeLinearOp_Def.hpp
  Piro_ObserverBase.hpp
  Piro_ROL_ObserverBase.hpp
  )
#
# B.2) epetra files
#
# Optional NOX Dependency (which includes LOCA)
IF (Piro_ENABLE_NOX)
  APPEND_SET(HEADERS
    Piro_SteadyStateSolver.hpp
    Piro_SteadyStateSolver_Def.hpp
    Piro_ProductModelEval.hpp
    Piro_NOXSolver.hpp
    Piro_NOXSolver_Def.hpp
    Piro_LOCASolver.hpp
    Piro_LOCASolver_Def.hpp
    Piro_ObserverToLOCASaveDataStrategyAdapter.hpp
    Piro_TrapezoidRuleSolver.hpp
    Piro_TrapezoidRuleSolver_Def.hpp
    Piro_TransientDecorator.hpp
  )
IF(PIRO_HAVE_EPETRA_STACK)
  APPEND_SET(HEADERS
    Piro_Epetra_NOXSolver.hpp
    Piro_Epetra_LOCASolver.hpp
    Piro_Epetra_LOCAAdaptiveSolver.hpp
    Piro_Epetra_AdaptiveSolutionManager.hpp
    Piro_Epetra_VelocityVerletSolver.hpp
    Piro_Epetra_NewmarkSolver.hpp
    Piro_Epetra_TrapezoidRuleSolver.hpp
    Piro_Epetra_SensitivityOperator.hpp
  )
ENDIF()
  APPEND_SET(SOURCES
    Piro_NOXSolver.cpp
    Piro_LOCASolver.cpp
    Piro_ObserverToLOCASaveDataStrategyAdapter.cpp
    )
IF(PIRO_HAVE_EPETRA_STACK)
  APPEND_SET(SOURCES
    Piro_Epetra_NOXSolver.cpp
    Piro_Epetra_LOCASolver.cpp
    Piro_Epetra_LOCAAdaptiveSolver.cpp
    Piro_Epetra_AdaptiveSolutionManager.cpp
    Piro_Epetra_VelocityVerletSolver.cpp
    Piro_Epetra_NewmarkSolver.cpp
    Piro_Epetra_TrapezoidRuleSolver.cpp
    Piro_Epetra_SensitivityOperator.cpp
    )
ENDIF()
ENDIF()
IF (Piro_ENABLE_Tempus)
  APPEND_SET(HEADERS
    Piro_TempusSolver.hpp
    Piro_TempusSolver_Def.hpp 
    Piro_TempusIntegrator.hpp
    Piro_TempusIntegrator_Def.hpp 
    Piro_TransientSolver.hpp
    Piro_TransientSolver_Def.hpp 
    Piro_TempusStepperFactory.hpp
    Piro_TempusStepControlFactory.hpp
    Piro_ObserverToTempusIntegrationObserverAdapter.hpp 
    Piro_ObserverToTempusIntegrationObserverAdapter_Def.hpp) 
ENDIF()

IF (Piro_ENABLE_ROL)
  APPEND_SET(HEADERS
    Piro_ThyraProductME_Objective_SimOpt.hpp
    Piro_ThyraProductME_Constraint_SimOpt.hpp
    Piro_CustomLBFGSSecant.hpp)

  # Add rol/adapters/thyra directories if Rol is enabled
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../rol/adapters/thyra/src/function)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../rol/adapters/thyra/src/vector)
ENDIF()
#

#
# B.2) Analysis Tools
#
APPEND_SET(HEADERS
  Piro_PerformSolve.hpp
  Piro_PerformSolve_Def.hpp
  )
APPEND_SET(SOURCES Piro_PerformSolve.cpp)

IF(PIRO_HAVE_EPETRA_STACK)
APPEND_SET(HEADERS Piro_Epetra_PerformSolve.hpp)
APPEND_SET(SOURCES Piro_Epetra_PerformSolve.cpp)
APPEND_SET(HEADERS Piro_Epetra_PerformAnalysis.hpp)
APPEND_SET(SOURCES Piro_Epetra_PerformAnalysis.cpp)
ENDIF()

APPEND_SET(HEADERS Piro_PerformAnalysis.hpp)
APPEND_SET(SOURCES Piro_PerformAnalysis.cpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  piro
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
