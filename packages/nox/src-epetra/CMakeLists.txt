
INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Base library
APPEND_SET(HEADERS
  NOX_Epetra.H
  NOX_Epetra_Vector.H
  NOX_Epetra_VectorSpace.H
  NOX_Epetra_VectorSpace_L2.H
  NOX_Epetra_VectorSpace_ScaledL2.H
  NOX_Epetra_MultiVector.H
  NOX_Epetra_Group.H
  NOX_Epetra_Interface_Required.H
  NOX_Epetra_Interface_Jacobian.H
  NOX_Epetra_Interface_Preconditioner.H
  NOX_Epetra_LinearSystem.H
  NOX_Epetra_LinearSystem_AztecOO.H
  NOX_Epetra_MatrixFree.H
  NOX_Epetra_FiniteDifference.H
  NOX_Epetra_FiniteDifferenceColoring.H
  NOX_Epetra_Scaling.H
  NOX_Epetra_BroydenOperator.H
  )

APPEND_SET(SOURCES
  NOX_Epetra_Vector.C
  NOX_Epetra_VectorSpace_L2.C
  NOX_Epetra_VectorSpace_ScaledL2.C
  NOX_Epetra_MultiVector.C
  NOX_Epetra_Group.C
  NOX_Epetra_LinearSystem_AztecOO.C
  NOX_Epetra_MatrixFree.C
  NOX_Epetra_FiniteDifference.C
  NOX_Epetra_FiniteDifferenceColoring.C
  NOX_Epetra_Scaling.C
  NOX_Epetra_BroydenOperator.C
  )

# Model Evaluator
IF (NOX_ENABLE_EpetraExt)
  APPEND_SET(HEADERS
    NOX_Epetra_ModelEvaluatorInterface.H
    )
  APPEND_SET(SOURCES
    NOX_Epetra_ModelEvaluatorInterface.C
    )
ENDIF()
     
# Amesos
IF (NOX_ENABLE_Amesos)
  APPEND_SET(HEADERS
    NOX_Epetra_LinearSystem_Amesos.H
    )
  APPEND_SET(SOURCES
    NOX_Epetra_LinearSystem_Amesos.C
    )
ENDIF()


IF (NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)
  PACKAGE_ADD_LIBRARY(
    noxepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS nox
    )
ENDIF()