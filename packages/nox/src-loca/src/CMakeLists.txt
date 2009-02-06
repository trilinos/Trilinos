
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(NOINSTALLHEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Base headers
APPEND_SET(HEADERS
  LOCA.H
  ) 

# Abstract
APPEND_SET(HEADERS
  LOCA_Abstract_Group.H
  LOCA_Abstract_TransposeSolveGroup.H
  LOCA_Abstract_Iterator.H
  LOCA_Abstract_Factory.H
  )

APPEND_SET(SOURCES
  LOCA_Abstract_Group.C
  LOCA_Abstract_Iterator.C
  LOCA_Abstract_Factory.C
  )

# Extended
APPEND_SET(HEADERS
  LOCA_Extended_Vector.H
  LOCA_Extended_MultiVector.H
  LOCA_Extended_MultiAbstractGroup.H
)
APPEND_SET(SOURCES
  LOCA_Extended_Vector.C
  LOCA_Extended_MultiVector.C
  LOCA_Extended_MultiAbstractGroup.C
)

# Bordered solver
APPEND_SET(HEADERS
  LOCA_BorderedSolver_AbstractOperator.H
  LOCA_BorderedSolver_JacobianOperator.H
  LOCA_BorderedSolver_ComplexOperator.H
  LOCA_BorderedSolver_BorderedOperator.H
  LOCA_BorderedSolver_AbstractStrategy.H
  LOCA_BorderedSolver_LowerTriangularBlockElimination.H
  LOCA_BorderedSolver_UpperTriangularBlockElimination.H
  LOCA_BorderedSolver_Bordering.H
  LOCA_BorderedSolver_HouseholderQR.H
  LOCA_BorderedSystem_AbstractGroup.H
  LOCA_BorderedSolver_Nested.H
  LOCA_BorderedSolver_Factory.H
  )
APPEND_SET(SOURCES
  LOCA_BorderedSolver_JacobianOperator.C
  LOCA_BorderedSolver_ComplexOperator.C
  LOCA_BorderedSolver_LowerTriangularBlockElimination.C
  LOCA_BorderedSolver_UpperTriangularBlockElimination.C
  LOCA_BorderedSolver_AbstractStrategy.C
  LOCA_BorderedSolver_Bordering.C
  LOCA_BorderedSolver_HouseholderQR.C
  LOCA_BorderedSolver_Nested.C
  LOCA_BorderedSolver_Factory.C
  )

# Multicontinuation
APPEND_SET(HEADERS
  LOCA_MultiContinuation_ExtendedVector.H
  LOCA_MultiContinuation_ExtendedMultiVector.H
  LOCA_MultiContinuation_AbstractGroup.H
  LOCA_MultiContinuation_FiniteDifferenceGroup.H
  LOCA_MultiContinuation_ConstraintInterface.H
  LOCA_MultiContinuation_ConstraintInterfaceMVDX.H
  LOCA_MultiContinuation_MultiVecConstraint.H
  LOCA_MultiContinuation_CompositeConstraint.H
  LOCA_MultiContinuation_CompositeConstraintMVDX.H
  LOCA_MultiContinuation_ConstrainedGroup.H
  LOCA_MultiContinuation_AbstractStrategy.H
  LOCA_MultiContinuation_ExtendedGroup.H
  LOCA_MultiContinuation_NaturalGroup.H
  LOCA_MultiContinuation_NaturalConstraint.H
  LOCA_MultiContinuation_ArcLengthGroup.H
  LOCA_MultiContinuation_ArcLengthConstraint.H
  LOCA_MultiContinuation_Factory.H 
  )
APPEND_SET(SOURCES
  LOCA_MultiContinuation_ExtendedVector.C
  LOCA_MultiContinuation_ExtendedMultiVector.C
  LOCA_MultiContinuation_AbstractGroup.C
  LOCA_MultiContinuation_FiniteDifferenceGroup.C
  LOCA_MultiContinuation_ConstraintInterfaceMVDX.C
  LOCA_MultiContinuation_MultiVecConstraint.C
  LOCA_MultiContinuation_CompositeConstraint.C
  LOCA_MultiContinuation_CompositeConstraintMVDX.C
  LOCA_MultiContinuation_ConstrainedGroup.C
  LOCA_MultiContinuation_ExtendedGroup.C
  LOCA_MultiContinuation_NaturalGroup.C
  LOCA_MultiContinuation_NaturalConstraint.C
  LOCA_MultiContinuation_ArcLengthGroup.C
  LOCA_MultiContinuation_ArcLengthConstraint.C
  LOCA_MultiContinuation_Factory.C
  )

# Bifurcation
APPEND_SET(HEADERS
  LOCA_TimeDependent_AbstractGroup.H
  LOCA_TurningPoint_MooreSpence_ExtendedVector.H
  LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.H
  LOCA_TurningPoint_MooreSpence_AbstractGroup.H
  LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H
  LOCA_TurningPoint_MooreSpence_ExtendedGroup.H
  LOCA_TurningPoint_MooreSpence_SolverStrategy.H
  LOCA_TurningPoint_MooreSpence_SalingerBordering.H
  LOCA_TurningPoint_MooreSpence_PhippsBordering.H
  LOCA_TurningPoint_MooreSpence_SolverFactory.H
  LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H
  LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H
  LOCA_TurningPoint_MinimallyAugmented_Constraint.H
  LOCA_TurningPoint_MinimallyAugmented_ModifiedConstraint.H
  LOCA_TurningPoint_MinimallyAugmented_ExtendedGroup.H
  LOCA_Pitchfork_MooreSpence_ExtendedVector.H
  LOCA_Pitchfork_MooreSpence_ExtendedMultiVector.H
  LOCA_Pitchfork_MooreSpence_AbstractGroup.H
  LOCA_Pitchfork_MooreSpence_ExtendedGroup.H
  LOCA_Pitchfork_MooreSpence_SolverStrategy.H
  LOCA_Pitchfork_MooreSpence_SalingerBordering.H
  LOCA_Pitchfork_MooreSpence_PhippsBordering.H
  LOCA_Pitchfork_MooreSpence_SolverFactory.H
  LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H
  LOCA_Pitchfork_MinimallyAugmented_Constraint.H
  LOCA_Pitchfork_MinimallyAugmented_ExtendedGroup.H
  LOCA_Hopf_MooreSpence_ExtendedVector.H
  LOCA_Hopf_MooreSpence_ExtendedMultiVector.H
  LOCA_Hopf_MooreSpence_AbstractGroup.H
  LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H
  LOCA_Hopf_MooreSpence_ExtendedGroup.H
  LOCA_Hopf_MooreSpence_SolverStrategy.H
  LOCA_Hopf_MooreSpence_SalingerBordering.H
  LOCA_Hopf_MooreSpence_SolverFactory.H
  LOCA_Hopf_ComplexVector.H
  LOCA_Hopf_ComplexMultiVector.H
  LOCA_Hopf_MinimallyAugmented_AbstractGroup.H
  LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H
  LOCA_Hopf_MinimallyAugmented_Constraint.H
  LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H
  LOCA_Bifurcation_Factory.H
  )
APPEND_SET(SOURCES
  LOCA_TurningPoint_MooreSpence_ExtendedVector.C
  LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.C
  LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.C
  LOCA_TurningPoint_MooreSpence_ExtendedGroup.C
  LOCA_TurningPoint_MooreSpence_SalingerBordering.C
  LOCA_TurningPoint_MooreSpence_PhippsBordering.C
  LOCA_TurningPoint_MooreSpence_SolverFactory.C
  LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.C
  LOCA_TurningPoint_MinimallyAugmented_Constraint.C
  LOCA_TurningPoint_MinimallyAugmented_ModifiedConstraint.C
  LOCA_TurningPoint_MinimallyAugmented_ExtendedGroup.C
  LOCA_Pitchfork_MooreSpence_ExtendedVector.C
  LOCA_Pitchfork_MooreSpence_ExtendedMultiVector.C
  LOCA_Pitchfork_MooreSpence_ExtendedGroup.C
  LOCA_Pitchfork_MooreSpence_SalingerBordering.C
  LOCA_Pitchfork_MooreSpence_PhippsBordering.C
  LOCA_Pitchfork_MooreSpence_SolverFactory.C
  LOCA_Pitchfork_MinimallyAugmented_Constraint.C
  LOCA_Pitchfork_MinimallyAugmented_ExtendedGroup.C
  LOCA_Hopf_MooreSpence_ExtendedVector.C
  LOCA_Hopf_MooreSpence_ExtendedMultiVector.C
  LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.C
  LOCA_Hopf_MooreSpence_ExtendedGroup.C
  LOCA_Hopf_MooreSpence_SalingerBordering.C
  LOCA_Hopf_MooreSpence_SolverFactory.C
  LOCA_Hopf_ComplexVector.C
  LOCA_Hopf_ComplexMultiVector.C
  LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.C
  LOCA_Hopf_MinimallyAugmented_Constraint.C
  LOCA_Hopf_MinimallyAugmented_ExtendedGroup.C
  LOCA_Bifurcation_Factory.C
  )

# Homotopy
APPEND_SET(HEADERS
  LOCA_Homotopy_AbstractGroup.H
        LOCA_Homotopy_Group.H
  LOCA_Homotopy_DeflatedGroup.H
  )
APPEND_SET(SOURCES
  LOCA_Homotopy_Group.C
  LOCA_Homotopy_DeflatedGroup.C
  )

# Utils
APPEND_SET(HEADERS
  LOCA_DerivUtils.H
  LOCA_ErrorCheck.H
  )
APPEND_SET(SOURCES
  LOCA_DerivUtils.C
  LOCA_ErrorCheck.C
  )

# Parameter
APPEND_SET(HEADERS
  LOCA_Parameter_Vector.H
  LOCA_Parameter_Library.H
  LOCA_Parameter_LibraryT.H
  LOCA_Parameter_Entry.H
  LOCA_Parameter_SublistParser.H
  )
APPEND_SET(SOURCES
  LOCA_Parameter_Vector.C
  LOCA_Parameter_Library.C
  LOCA_Parameter_SublistParser.C
  )

# Predictor
APPEND_SET(HEADERS
  LOCA_MultiPredictor_AbstractStrategy.H
  LOCA_MultiPredictor_Constant.H
  LOCA_MultiPredictor_Tangent.H
  LOCA_MultiPredictor_Secant.H
  LOCA_MultiPredictor_Random.H
  LOCA_MultiPredictor_Restart.H
  LOCA_MultiPredictor_Factory.H
  )
APPEND_SET(SOURCES
  LOCA_MultiPredictor_AbstractStrategy.C
  LOCA_MultiPredictor_Constant.C
  LOCA_MultiPredictor_Tangent.C
  LOCA_MultiPredictor_Secant.C
  LOCA_MultiPredictor_Random.C
  LOCA_MultiPredictor_Restart.C
  LOCA_MultiPredictor_Factory.C
  )

# Stepsize
APPEND_SET(HEADERS
  LOCA_StepSize_AbstractStrategy.H
  LOCA_StepSize_Constant.H
  LOCA_StepSize_Adaptive.H
  LOCA_StepSize_Factory.H
  )
APPEND_SET(SOURCES
  LOCA_StepSize_Constant.C
  LOCA_StepSize_Adaptive.C
  LOCA_StepSize_Factory.C
  )

# Stepper
APPEND_SET(HEADERS
  LOCA_Stepper.H
  )
APPEND_SET(SOURCES
  LOCA_Stepper.C
  )

# Status test
APPEND_SET(HEADERS
  LOCA_Solver_Wrapper.H
  LOCA_StatusTest_Wrapper.H
  )
APPEND_SET(SOURCES
  LOCA_Solver_Wrapper.C
  LOCA_StatusTest_Wrapper.C
  )

# Factory
APPEND_SET(HEADERS
  LOCA_Factory.H
  LOCA_GlobalData.H
  )
APPEND_SET(SOURCES
  LOCA_Factory.C
  LOCA_GlobalData.C
  )

# Eigensolver
APPEND_SET(HEADERS
  LOCA_Eigensolver_AbstractStrategy.H
  LOCA_Eigensolver_DefaultStrategy.H
  LOCA_Eigensolver_Factory.H
  LOCA_EigenvalueSort_Strategies.H
  LOCA_EigenvalueSort_Factory.H
  LOCA_SaveEigenData_AbstractStrategy.H
  LOCA_SaveEigenData_DefaultStrategy.H
  LOCA_SaveEigenData_Factory.H
  )
APPEND_SET(SOURCES
  LOCA_Eigensolver_DefaultStrategy.C
  LOCA_Eigensolver_Factory.C
  LOCA_EigenvalueSort_Strategies.C
  LOCA_EigenvalueSort_Factory.C
  LOCA_SaveEigenData_DefaultStrategy.C
  LOCA_SaveEigenData_Factory.C
  )

# Anasazi
APPEND_SET(HEADERS
  LOCA_AnasaziOperator_AbstractStrategy.H
  LOCA_AnasaziOperator_JacobianInverse.H
  LOCA_AnasaziOperator_ShiftInvert.H
  LOCA_AnasaziOperator_Cayley.H
  LOCA_AnasaziOperator_Factory.H
  )
APPEND_SET(SOURCES
  LOCA_AnasaziOperator_JacobianInverse.C
  LOCA_AnasaziOperator_ShiftInvert.C
  LOCA_AnasaziOperator_Cayley.C
  LOCA_AnasaziOperator_Factory.C
  )
IF (NOX_ENABLE_Anasazi)
  APPEND_SET(HEADERS
    Anasazi_LOCA_MultiVecTraits.H
    Anasazi_LOCA_OperatorTraits.H
    Anasazi_LOCA_Sort.H
    LOCA_Eigensolver_AnasaziStrategy.H
    )
  APPEND_SET(SOURCES
    Anasazi_LOCA_Sort.C
    LOCA_Eigensolver_AnasaziStrategy.C
    )
ENDIF()

#
# C) Define the targets for package's library(s)
#

IF (NOX_ENABLE_LOCA)
  PACKAGE_ADD_LIBRARY(
    loca
    HEADERS ${HEADERS}
    NOINSTALLHEADERS ${NOINSTALLHEADERS}
    SOURCES ${SOURCES}
    )
ENDIF()
