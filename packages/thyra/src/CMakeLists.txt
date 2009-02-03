
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  ${PACKAGE_NAME}_ConfigDefs.hpp 
  ${PACKAGE_NAME}_Version.hpp
  )

APPEND_SET(SOURCES
  ${PACKAGE_NAME}_Version.cpp
  )

#
# src/interfaces/operator_vector/ana/fundamental
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_vector/ana/fundamental)

APPEND_SET(HEADERS
  interfaces/operator_vector/ana/fundamental/Thyra_LinearOpBaseDecl.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_LinearOpBase.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_MultiVectorBaseDecl.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_MultiVectorBase.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_OperatorVectorTypes.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_VectorBaseDecl.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_VectorBase.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_VectorSpaceBaseDecl.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_VectorSpaceBase.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_VectorSpaceFactoryBaseDecl.hpp
  interfaces/operator_vector/ana/fundamental/Thyra_VectorSpaceFactoryBase.hpp
  )

#
# src/interfaces/operator_vector/ana/extended
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_vector/ana/extended)

APPEND_SET(HEADERS
  interfaces/operator_vector/ana/extended/Thyra_AddedLinearOpBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_BlockedLinearOpBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_DiagonalLinearOpBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_IdentityLinearOpBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_MultipliedLinearOpBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_MultiVectorFileIOBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_MultiVectorRandomizerBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_PhysicallyBlockedLinearOpBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_PolynomialVectorTraits.hpp
  interfaces/operator_vector/ana/extended/Thyra_ProductMultiVectorBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_ProductVectorBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_ProductVectorSpaceBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_ScaledAdjointLinearOpBaseDecl.hpp
  interfaces/operator_vector/ana/extended/Thyra_ScaledAdjointLinearOpBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_SpmdVectorSpaceBaseDecl.hpp
  interfaces/operator_vector/ana/extended/Thyra_SpmdVectorSpaceBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_VectorSpaceConverterBase.hpp
  interfaces/operator_vector/ana/extended/Thyra_ZeroLinearOpBase.hpp
  )

#
# src/support/operator_vector/client_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/client_support)

APPEND_SET(HEADERS
  support/operator_vector/client_support/Thyra_AssertOp.hpp
  support/operator_vector/client_support/Thyra_DefaultAddedLinearOp.hpp
  support/operator_vector/client_support/Thyra_DefaultAddedLinearOpDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultBlockedLinearOp.hpp
  support/operator_vector/client_support/Thyra_DefaultBlockedLinearOpDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultColumnwiseMultiVector.hpp
  support/operator_vector/client_support/Thyra_DefaultColumnwiseMultiVectorDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultDiagonalLinearOp.hpp
  support/operator_vector/client_support/Thyra_DefaultDiagonalLinearOpDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultIdentityLinearOp.hpp
  support/operator_vector/client_support/Thyra_DefaultIdentityLinearOpDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultMultiVectorProductVector.hpp
  support/operator_vector/client_support/Thyra_DefaultMultiVectorProductVectorDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultMultiVectorProductVectorSpace.hpp
  support/operator_vector/client_support/Thyra_DefaultMultiVectorProductVectorSpaceDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultMultipliedLinearOp.hpp
  support/operator_vector/client_support/Thyra_DefaultMultipliedLinearOpDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultProductMultiVector.hpp
  support/operator_vector/client_support/Thyra_DefaultProductMultiVectorDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultProductVector.hpp
  support/operator_vector/client_support/Thyra_DefaultProductVectorDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultProductVectorSpace.hpp
  support/operator_vector/client_support/Thyra_DefaultProductVectorSpaceDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultScaledAdjointLinearOp.hpp
  support/operator_vector/client_support/Thyra_DefaultScaledAdjointLinearOpDecl.hpp
  support/operator_vector/client_support/Thyra_DefaultZeroLinearOp.hpp
  support/operator_vector/client_support/Thyra_DefaultZeroLinearOpDecl.hpp
  support/operator_vector/client_support/Thyra_DetachedMultiVectorView.hpp
  support/operator_vector/client_support/Thyra_DetachedVectorView.hpp
  support/operator_vector/client_support/Thyra_LinearCombinationDecl.hpp
  support/operator_vector/client_support/Thyra_LinearCombinationImpl.hpp
  support/operator_vector/client_support/Thyra_LinearCombinationTester.hpp
  support/operator_vector/client_support/Thyra_LinearOpDefaultBase.hpp
  support/operator_vector/client_support/Thyra_LinearOpDefaultBaseDecl.hpp
  support/operator_vector/client_support/Thyra_LinearOpTester.hpp
  support/operator_vector/client_support/Thyra_LinearOpTesterDecl.hpp
  support/operator_vector/client_support/Thyra_LinearOperatorDecl.hpp
  support/operator_vector/client_support/Thyra_LinearOperatorImpl.hpp
  support/operator_vector/client_support/Thyra_ListedMultiVectorRandomizer.hpp
  support/operator_vector/client_support/Thyra_MultiVectorAllocator.hpp
  support/operator_vector/client_support/Thyra_MultiVectorDefaultBase.hpp
  support/operator_vector/client_support/Thyra_MultiVectorDefaultBaseDecl.hpp
  support/operator_vector/client_support/Thyra_MultiVectorStdOps.hpp
  support/operator_vector/client_support/Thyra_MultiVectorStdOpsDecl.hpp
  support/operator_vector/client_support/Thyra_MultiVectorStdOpsTester.hpp
  support/operator_vector/client_support/Thyra_MultiVectorStdOpsTesterDecl.hpp
  support/operator_vector/client_support/Thyra_MultiVectorTester.hpp
  support/operator_vector/client_support/Thyra_MultiVectorTesterDecl.hpp
  support/operator_vector/client_support/Thyra_OperatorVectorClientSupport.hpp
  support/operator_vector/client_support/Thyra_ParameterDrivenMultiVectorInput.hpp
  support/operator_vector/client_support/Thyra_SingleRhsLinearOpBase.hpp
  support/operator_vector/client_support/Thyra_SingleRhsLinearOpBaseDecl.hpp
  support/operator_vector/client_support/Thyra_SingleScalarLinearOpBase.hpp
  support/operator_vector/client_support/Thyra_SingleScalarLinearOpBaseDecl.hpp
  support/operator_vector/client_support/Thyra_TestSpecifier.hpp
  support/operator_vector/client_support/Thyra_TesterBase.hpp
  support/operator_vector/client_support/Thyra_TestingTools.hpp
  support/operator_vector/client_support/Thyra_TestingToolsDecl.hpp
  support/operator_vector/client_support/Thyra_UniversalMultiVectorRandomizer.hpp
  support/operator_vector/client_support/Thyra_VecOpMacros.hpp
  support/operator_vector/client_support/Thyra_VectorDecl.hpp
  support/operator_vector/client_support/Thyra_VectorDefaultBase.hpp
  support/operator_vector/client_support/Thyra_VectorDefaultBaseDecl.hpp
  support/operator_vector/client_support/Thyra_VectorHandleOpsDecl.hpp
  support/operator_vector/client_support/Thyra_VectorHandleOpsImpl.hpp
  support/operator_vector/client_support/Thyra_VectorImpl.hpp
  support/operator_vector/client_support/Thyra_VectorOpTester.hpp
  support/operator_vector/client_support/Thyra_VectorSpaceDecl.hpp
  support/operator_vector/client_support/Thyra_VectorSpaceDefaultBase.hpp
  support/operator_vector/client_support/Thyra_VectorSpaceDefaultBaseDecl.hpp
  support/operator_vector/client_support/Thyra_VectorSpaceImpl.hpp
  support/operator_vector/client_support/Thyra_VectorSpaceTester.hpp
  support/operator_vector/client_support/Thyra_VectorSpaceTesterDecl.hpp
  support/operator_vector/client_support/Thyra_VectorStdOps.hpp
  support/operator_vector/client_support/Thyra_VectorStdOpsDecl.hpp
  support/operator_vector/client_support/Thyra_VectorStdOpsTester.hpp
  support/operator_vector/client_support/Thyra_VectorStdOpsTesterDecl.hpp
  support/operator_vector/client_support/Thyra_VectorTester.hpp
  support/operator_vector/client_support/Thyra_VectorTesterDecl.hpp
  support/operator_vector/client_support/Thyra_describeLinearOp.hpp
  support/operator_vector/client_support/Thyra_describeLinearOpDecl.hpp
  )

APPEND_SET(SOURCES
  support/operator_vector/client_support/Thyra_TestingTools.cpp
  )

#
# src/support/operator_vector/adapter_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/adapter_support)

APPEND_SET(HEADERS
  support/operator_vector/adapter_support/Thyra_DefaultClusteredSpmdProductVector.hpp
  support/operator_vector/adapter_support/Thyra_DefaultClusteredSpmdProductVectorDecl.hpp
  support/operator_vector/adapter_support/Thyra_DefaultClusteredSpmdProductVectorSpace.hpp
  support/operator_vector/adapter_support/Thyra_DefaultClusteredSpmdProductVectorSpaceDecl.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdMultiVector.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdMultiVectorDecl.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdMultiVectorFileIO.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdVector.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorDecl.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorSpace.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorSpaceDecl.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorSpaceFactory.hpp
  support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorSpaceFactoryDecl.hpp
  support/operator_vector/adapter_support/Thyra_EuclideanLinearOpBase.hpp
  support/operator_vector/adapter_support/Thyra_EuclideanLinearOpBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_EuclideanScalarProd.hpp
  support/operator_vector/adapter_support/Thyra_EuclideanScalarProdDecl.hpp
  support/operator_vector/adapter_support/Thyra_LinearOpScalarProd.hpp
  support/operator_vector/adapter_support/Thyra_LinearOpScalarProdDecl.hpp
  support/operator_vector/adapter_support/Thyra_OperatorVectorAdapterSupportTypes.hpp
  support/operator_vector/adapter_support/Thyra_ScalarProdBase.hpp
  support/operator_vector/adapter_support/Thyra_ScalarProdBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_ScalarProdVectorSpaceBase.hpp
  support/operator_vector/adapter_support/Thyra_ScalarProdVectorSpaceBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_SingleRhsEuclideanLinearOpBase.hpp
  support/operator_vector/adapter_support/Thyra_SingleRhsEuclideanLinearOpBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_SingleScalarEuclideanLinearOpBase.hpp
  support/operator_vector/adapter_support/Thyra_SingleScalarEuclideanLinearOpBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_SpmdLinearOpBase.hpp
  support/operator_vector/adapter_support/Thyra_SpmdLinearOpBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_SpmdMultiVectorBase.hpp
  support/operator_vector/adapter_support/Thyra_SpmdMultiVectorBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_SpmdMultiVectorSerializer.hpp
  support/operator_vector/adapter_support/Thyra_SpmdMultiVectorSerializerDecl.hpp
  support/operator_vector/adapter_support/Thyra_SpmdVectorBase.hpp
  support/operator_vector/adapter_support/Thyra_SpmdVectorBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_SpmdVectorSpaceDefaultBase.hpp
  support/operator_vector/adapter_support/Thyra_SpmdVectorSpaceDefaultBaseDecl.hpp
  support/operator_vector/adapter_support/Thyra_SpmdVectorSpaceUtilities.hpp
  support/operator_vector/adapter_support/Thyra_apply_op_helper.hpp
  support/operator_vector/adapter_support/Thyra_apply_op_helper_decl.hpp
  )

APPEND_SET(SOURCES
  support/operator_vector/adapter_support/Thyra_SpmdVectorSpaceUtilities.cpp
  )

#
# src/support/operator_vector/adapter_support/serial
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/adapter_support/serial)

APPEND_SET(HEADERS
  support/operator_vector/adapter_support/serial/Thyra_DefaultSerialVectorSpaceConverterDecl.hpp
  support/operator_vector/adapter_support/serial/Thyra_DefaultSerialVectorSpaceConverter.hpp
  support/operator_vector/adapter_support/serial/Thyra_SerialLinearOpBaseDecl.hpp
  support/operator_vector/adapter_support/serial/Thyra_SerialLinearOpBase.hpp
  support/operator_vector/adapter_support/serial/Thyra_SerialVectorSpaceConverterBaseDecl.hpp
  support/operator_vector/adapter_support/serial/Thyra_SerialVectorSpaceConverterBase.hpp
  )

#
# src/interfaces/operator_solve/ana/fundamental
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_solve/ana/fundamental)

APPEND_SET(HEADERS
  interfaces/operator_solve/ana/fundamental/Thyra_LinearOpSourceBase.hpp
  interfaces/operator_solve/ana/fundamental/Thyra_LinearOpWithSolveBase.hpp
  interfaces/operator_solve/ana/fundamental/Thyra_LinearOpWithSolveBaseDecl.hpp
  interfaces/operator_solve/ana/fundamental/Thyra_LinearOpWithSolveFactoryBase.hpp
  interfaces/operator_solve/ana/fundamental/Thyra_PreconditionerBase.hpp
  interfaces/operator_solve/ana/fundamental/Thyra_PreconditionerFactoryBase.hpp
  interfaces/operator_solve/ana/fundamental/Thyra_SolveSupportTypes.hpp
  )

#
# src/interfaces/operator_solve/ana/extended
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_solve/ana/extended)

APPEND_SET(HEADERS
  interfaces/operator_solve/ana/extended/Thyra_BlockedLinearOpWithSolveBase.hpp
  interfaces/operator_solve/ana/extended/Thyra_InverseLinearOpBase.hpp
  interfaces/operator_solve/ana/extended/Thyra_LinearSolverBuilderBase.hpp
  interfaces/operator_solve/ana/extended/Thyra_PhysicallyBlockedLinearOpWithSolveBase.hpp
  )

#
# src/support/operator_solve/client_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/operator_solve/client_support)

APPEND_SET(HEADERS
  support/operator_solve/client_support/Thyra_DefaultAdjointLinearOpWithSolve.hpp
  support/operator_solve/client_support/Thyra_DefaultAdjointLinearOpWithSolveDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultBlockedTriangularLinearOpWithSolve.hpp
  support/operator_solve/client_support/Thyra_DefaultBlockedTriangularLinearOpWithSolveDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultBlockedTriangularLinearOpWithSolveFactory.hpp
  support/operator_solve/client_support/Thyra_DefaultBlockedTriangularLinearOpWithSolveFactoryDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultDiagonalLinearOpWithSolve.hpp
  support/operator_solve/client_support/Thyra_DefaultDiagonalLinearOpWithSolveDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultInverseLinearOp.hpp
  support/operator_solve/client_support/Thyra_DefaultInverseLinearOpDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultLinearOpSource.hpp
  support/operator_solve/client_support/Thyra_DefaultLinearOpSourceDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultMultiVectorLinearOpWithSolve.hpp
  support/operator_solve/client_support/Thyra_DefaultMultiVectorLinearOpWithSolveDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultPreconditioner.hpp
  support/operator_solve/client_support/Thyra_DefaultPreconditionerDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultSerialDenseLinearOpWithSolve.hpp
  support/operator_solve/client_support/Thyra_DefaultSerialDenseLinearOpWithSolveDecl.hpp
  support/operator_solve/client_support/Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp
  support/operator_solve/client_support/Thyra_DefaultSerialDenseLinearOpWithSolveFactoryDecl.hpp
  support/operator_solve/client_support/Thyra_DelayedLinearOpWithSolve.hpp
  support/operator_solve/client_support/Thyra_DelayedLinearOpWithSolveDecl.hpp
  support/operator_solve/client_support/Thyra_DelayedLinearOpWithSolveFactory.hpp
  support/operator_solve/client_support/Thyra_DelayedLinearOpWithSolveFactoryDecl.hpp
  support/operator_solve/client_support/Thyra_InverseLinearOperator.hpp
  support/operator_solve/client_support/Thyra_LinearOpWithSolveFactoryHelpers.hpp
  support/operator_solve/client_support/Thyra_LinearOpWithSolveHelpers.hpp
  support/operator_solve/client_support/Thyra_LinearOpWithSolveTester.hpp
  support/operator_solve/client_support/Thyra_LinearOpWithSolveTesterDecl.hpp
  support/operator_solve/client_support/Thyra_PreconditionerFactoryHelpers.hpp
  support/operator_solve/client_support/Thyra_SingleRhsLinearOpWithSolveBase.hpp
  support/operator_solve/client_support/Thyra_SingleRhsLinearOpWithSolveBaseDecl.hpp
  support/operator_solve/client_support/Thyra_SingleScalarLinearOpWithSolveBase.hpp
  support/operator_solve/client_support/Thyra_SingleScalarLinearOpWithSolveBaseDecl.hpp
  )

#
# src/interfaces/nonlinear/model_evaluator/ana/fundamental
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/nonlinear/model_evaluator/ana/fundamental)

APPEND_SET(HEADERS
  interfaces/nonlinear/model_evaluator/ana/fundamental/Thyra_ModelEvaluator.hpp
  interfaces/nonlinear/model_evaluator/ana/fundamental/Thyra_ModelEvaluatorBase.hpp
  interfaces/nonlinear/model_evaluator/ana/fundamental/Thyra_ModelEvaluatorBase_decl.hpp
  interfaces/nonlinear/model_evaluator/ana/fundamental/Thyra_ModelEvaluatorBase_def.hpp
  )

APPEND_SET(SOURCES
  interfaces/nonlinear/model_evaluator/ana/fundamental/Thyra_ModelEvaluatorBase.cpp
  )

#
# src/support/nonlinear/model_evaluator/client_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/nonlinear/model_evaluator/client_support)

APPEND_SET(HEADERS
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultLumpedParameterModelEvaluator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultEvaluationLoggerModelEvaluator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultFinalPointCaptureModelEvaluator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultFiniteDifferenceModelEvaluator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultInverseModelEvaluator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultLumpedParameterModelEvaluator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultModelEvaluatorWithSolveFactory.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultMultiPeriodModelEvaluator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_DirectionalFiniteDiffCalculator.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_ModelEvaluatorDefaultBase.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_ModelEvaluatorDelegatorBase.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_ModelEvaluatorHelpers.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_ResponseOnlyModelEvaluatorBase.hpp
  support/nonlinear/model_evaluator/client_support/Thyra_StateFuncModelEvaluatorBase.hpp
  )

#
# src/interfaces/nonlinear/solvers/ana/fundamental
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/nonlinear/solvers/ana/fundamental)

APPEND_SET(HEADERS
  interfaces/nonlinear/solvers/ana/fundamental/Thyra_NonlinearSolverBase.hpp
  )

#
# src/support/nonlinear/solvers/client_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/nonlinear/solvers/client_support)

APPEND_SET(HEADERS
  support/nonlinear/solvers/client_support/Thyra_DampenedNewtonNonlinearSolver.hpp
  support/nonlinear/solvers/client_support/Thyra_DefaultStateEliminationModelEvaluator.hpp
  support/nonlinear/solvers/client_support/Thyra_LinearNonlinearSolver.hpp
  )

#
# thyra/example/operator_vector
#

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/example/operator_vector)

APPEND_SET(HEADERS
  ${PACKAGE_SOURCE_DIR}/example/operator_vector/sillyModifiedGramSchmidt.hpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  thyra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
