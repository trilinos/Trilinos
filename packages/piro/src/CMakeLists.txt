
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
  Piro_MatrixFreeDecorator.hpp
  Piro_MatrixFreeDecorator_Def.hpp
  Piro_MatrixFreeLinearOp.hpp
  Piro_MatrixFreeLinearOp_Def.hpp
  Piro_ObserverBase.hpp
  Piro_ROL_ObserverBase.hpp
)

APPEND_SET(SOURCES
  Piro_ValidPiroParameters.cpp
  Piro_StratimikosUtils.cpp
  Piro_Version.cpp
)
#
# B.2) solvers
#
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
  APPEND_SET(SOURCES
    Piro_NOXSolver.cpp
    Piro_LOCASolver.cpp
    Piro_ObserverToLOCASaveDataStrategyAdapter.cpp
    )
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
  TRIBITS_INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../rol/adapters/thyra/src/function
    ${CMAKE_CURRENT_SOURCE_DIR}/../../rol/adapters/thyra/src/vector
  )

  IF (Piro_ENABLE_HdsaLib)
    APPEND_SET(HEADERS
      Piro_HDSA_MD_ROL_Data_Interface.hpp
      Piro_HDSA_MD_ROL_Elliptic_u_Prior_Interface.hpp
      Piro_HDSA_MD_ROL_Elliptic_z_Prior_Interface.hpp)
      
    TRIBITS_INCLUDE_DIRECTORIES(
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/interfaces/hyperparameter_interfaces
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/interfaces/problem_interfaces
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/interfaces/ouu_interfaces
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/interfaces/transient_prior_interfaces
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/interfaces/u_prior_interfaces
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/interfaces/z_prior_interfaces
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/analysis/prior
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/analysis/posterior
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/model_discrepancy/analysis/auxiliary
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/pseudo_time_continuation/interfaces
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/pseudo_time_continuation/analysis
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/base/auxiliary
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/base/linear_algebra
      ${HDSALIB_REPO_SOURCE_DIR}/src/core/base/vectors
      ${HDSALIB_REPO_SOURCE_DIR}/src/trilinos_adapter/adapters
      ${HDSALIB_REPO_SOURCE_DIR}/src/trilinos_adapter/core
      ${HDSALIB_REPO_SOURCE_DIR}/interfaces/rol
    )
  ENDIF()
ENDIF()
#

#
# B.2) drivers
#
APPEND_SET(HEADERS
  Piro_PerformSolve.hpp
  Piro_PerformSolve_Def.hpp
  Piro_PerformAnalysis.hpp
)
APPEND_SET(SOURCES 
  Piro_PerformSolve.cpp
  Piro_PerformAnalysis.cpp
)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  piro
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
