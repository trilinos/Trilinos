
INCLUDE(TribitsLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/continuation-manager)

# Base library
APPEND_SET(HEADERS
  LOCA_Epetra.H
  LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H
  LOCA_Epetra_TransposeLinearSystem_TransposePreconditioner.H
  LOCA_Epetra_IdentityOp.H
  LOCA_Epetra_LeftPreconditionedOp.H
  LOCA_Epetra_TransposeLinearSystem_LeftPreconditioning.H
  LOCA_Epetra_TransposeLinearSystem_Factory.H
  LOCA_Epetra_Interface_Required.H
  LOCA_Epetra_Interface_MassMatrix.H
  LOCA_Epetra_Interface_TimeDependent.H
  LOCA_Epetra_Interface_TimeDependentMatrixFree.H
  LOCA_Epetra_Interface_FreeEnergy.H
  LOCA_Epetra_Group.H
  LOCA_Epetra_CompactWYOp.H
  LOCA_Epetra_LowRankUpdateOp.H
  LOCA_Epetra_LowRankUpdateRowMatrix.H
  LOCA_BorderedSolver_EpetraHouseholder.H
  LOCA_Epetra_AugmentedOp.H
  LOCA_BorderedSolver_EpetraAugmented.H
  LOCA_Epetra_Factory.H
  LOCA_Epetra_ShiftInvertOperator.H
  LOCA_Epetra_AdaptiveSolutionManager.H
  )

APPEND_SET(SOURCES
  LOCA_Epetra_TransposeLinearSystem_TransposePreconditioner.C
  LOCA_Epetra_IdentityOp.C
  LOCA_Epetra_LeftPreconditionedOp.C
  LOCA_Epetra_TransposeLinearSystem_LeftPreconditioning.C
  LOCA_Epetra_TransposeLinearSystem_Factory.C
  LOCA_Epetra_Group.C
  LOCA_Epetra_CompactWYOp.C
  LOCA_Epetra_LowRankUpdateOp.C
  LOCA_Epetra_LowRankUpdateRowMatrix.C
  LOCA_BorderedSolver_EpetraHouseholder.C
  LOCA_Epetra_AugmentedOp.C
  LOCA_BorderedSolver_EpetraAugmented.C
  LOCA_Epetra_Factory.C
  LOCA_Epetra_ShiftInvertOperator.C
  LOCA_Epetra_AdaptiveSolutionManager.C
  )

# EpetraExt stuff
IF (NOX_ENABLE_EpetraExt)
  APPEND_SET(HEADERS
    LOCA_Epetra_Interface_xyzt.H
    LOCA_Epetra_Interface_MultiPoint.H
    LOCA_Epetra_xyztPrec.H
    LOCA_Epetra_AnasaziOperator_Floquet.H
    LOCA_Epetra_TransposeLinearSystem_ExplicitTranspose.H
    LOCA_Epetra_ModelEvaluatorInterface.H
    LOCA_Epetra_AdaptiveStepper.H
    continuation-manager/ContinuationManager.H
    continuation-manager/LOCAInterface.H
    continuation-manager/ProblemNOXPrototype.H
    continuation-manager/ProblemLOCAPrototype.H
    continuation-manager/IOContFileUtils.H
    )
  APPEND_SET(SOURCES
    LOCA_Epetra_Interface_xyzt.C
    LOCA_Epetra_Interface_MultiPoint.C
    LOCA_Epetra_xyztPrec.C
    LOCA_Epetra_AnasaziOperator_Floquet.C
    LOCA_Epetra_TransposeLinearSystem_ExplicitTranspose.C
    LOCA_Epetra_ModelEvaluatorInterface.C
    LOCA_Epetra_AdaptiveStepper.C
    continuation-manager/ContinuationManager.C
    continuation-manager/LOCAInterface.C
    continuation-manager/ProblemNOXPrototype.C
    continuation-manager/ProblemLOCAPrototype.C
    continuation-manager/IOContFileUtils.C
    )
ENDIF()

IF (NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)
  TRIBITS_ADD_LIBRARY(
    locaepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS loca noxepetra nox
    )
ENDIF()
