
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

ASSERT_DEFINED(Teuchos_ENABLE_Boost)
SET(HAVE_RYTHMOS_TEUCHOS_BOOST ${Teuchos_ENABLE_Boost})

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Rythmos_AdjointModelEvaluator.hpp
  Rythmos_BackwardEulerStepper.hpp
  Rythmos_BreakPointInformerBase.hpp
  Rythmos_CompositeIntegratorObserver.hpp
  Rythmos_ConfigDefs.h
  Rythmos_CubicSplineInterpolator.hpp
  Rythmos_DataStore.hpp
  Rythmos_DefaultIntegrator.hpp
  Rythmos_DiagonalImplicitRKModelEvaluator.hpp
  Rythmos_ErrorEstimateBase.hpp
  Rythmos_ErrorEstimatorBase.hpp
  Rythmos_ErrWtVecCalcAcceptingStepControlStrategyBase.hpp
  Rythmos_ErrWtVecCalcBase.hpp
  Rythmos_ExplicitRKStepper.hpp
  Rythmos_ExplicitTaylorPolynomialStepper.hpp
  Rythmos_ForwardEulerStepper.hpp
  Rythmos_ForwardResponseSensitivityComputer.hpp
  Rythmos_ForwardResponseSensitivityComputerObserver.hpp
  Rythmos_ForwardSensitivityIntegratorAsModelEvaluator.hpp
  Rythmos_ForwardSensitivityModelEvaluator.hpp
  Rythmos_ForwardSensitivityStepper.hpp
  Rythmos_HermiteInterpolator.hpp
  Rythmos_ImplicitBDFStepper.hpp
  Rythmos_ImplicitBDFStepperDecl.hpp
  Rythmos_ImplicitBDFStepperErrWtVecCalc.hpp
  Rythmos_ImplicitBDFStepperStepControl.hpp
  Rythmos_ImplicitBDFStepperStepControlDecl.hpp
  Rythmos_ImplicitRKModelEvaluator.hpp
  Rythmos_ImplicitRKStepper.hpp
  Rythmos_IntegrationControlStrategyBase.hpp
  Rythmos_IntegrationObserverBase.hpp
  Rythmos_IntegratorBase.hpp
  Rythmos_InterpolationBuffer.hpp
  Rythmos_InterpolationBufferAppenderBase.hpp
  Rythmos_InterpolationBufferBase.hpp
  Rythmos_InterpolationBufferHelpers.hpp
  Rythmos_InterpolatorBase.hpp
  Rythmos_InterpolatorBaseHelpers.hpp
  Rythmos_LinearInterpolator.hpp
  Rythmos_PointwiseInterpolationBufferAppender.hpp
  Rythmos_QuadratureBase.hpp
  Rythmos_RKButcherTableau.hpp
  Rythmos_ResponseAndFwdSensPoint.hpp
  Rythmos_SimpleIntegrationControlStrategy.hpp
  Rythmos_SingleResidualModelEvaluator.hpp
  Rythmos_SingleResidualModelEvaluatorBase.hpp
  Rythmos_SmartInterpolationBufferAppender.hpp
  Rythmos_SolverAcceptingStepperBase.hpp
  Rythmos_StateAndForwardSensitivityModelEvaluator.hpp
  Rythmos_StepControlInfo.hpp
  Rythmos_StepControlStrategyAcceptingStepperBase.hpp
  Rythmos_StepControlStrategyBase.hpp
  Rythmos_StepperAsModelEvaluator.hpp
  Rythmos_StepperBase.hpp
  Rythmos_StepperBuilder.hpp
  Rythmos_StepperHelpers.hpp
  Rythmos_StepperSupportTypes.hpp
  Rythmos_TimeDiscretizedBackwardEulerModelEvaluator.hpp
  Rythmos_TimeRange.hpp
  Rythmos_TimeStepNonlinearSolver.hpp
  Rythmos_Types.hpp
  Rythmos_extractStateAndSens.hpp
  Rythmos_Version.h
  )

APPEND_SET(SOURCES
  Rythmos_Version.cpp
  Rythmos_RKButcherTableau.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  rythmos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


#
# Add a test-only model library
#

SET(HEADERS)
SET(SOURCES)

INCLUDE_DIRECTORIES( 
  ${PACKAGE_SOURCE_DIR}/test/SinCos
  ${PACKAGE_SOURCE_DIR}/test/PolynomialModel
  )

APPEND_SET(HEADERS
  ${PACKAGE_SOURCE_DIR}/test/SinCos/SinCosModel.hpp
  ${PACKAGE_SOURCE_DIR}/test/PolynomialModel/PolynomialModel.hpp
  )

APPEND_SET(SOURCES
  ${PACKAGE_SOURCE_DIR}/test/SinCos/SinCosModel.cpp
  ${PACKAGE_SOURCE_DIR}/test/PolynomialModel/PolynomialModel.cpp
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_EpetraExt)

  INCLUDE_DIRECTORIES( 
    ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalTransient
    )
  
  APPEND_SET(HEADERS
    ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.hpp
    )
  
  APPEND_SET(SOURCES
    ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp
    )

ENDIF()

PACKAGE_ADD_LIBRARY(
  rythmos_test_models
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  DEPLIBS rythmos
  )




