
INCLUDE(TrilinosPackageLibraryHelpers)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS
  Thyra_DiagonalEpetraLinearOpWithSolveFactory.hpp
  Thyra_EpetraTypes.hpp
  Thyra_EpetraLinearOpBase.hpp
  Thyra_EpetraLinearOp.hpp
  Thyra_EpetraOperatorViewExtractorBase.hpp
  Thyra_EpetraOperatorViewExtractorStd.hpp
  Thyra_EpetraOperatorWrapper.hpp
  Thyra_EpetraThyraWrappers.hpp
  Thyra_get_Epetra_Operator.hpp
  )

SET(SOURCES
  Thyra_DiagonalEpetraLinearOpWithSolveFactory.cpp
  Thyra_EpetraLinearOp.cpp
  Thyra_EpetraOperatorViewExtractorStd.cpp
  Thyra_EpetraOperatorWrapper.cpp
  Thyra_EpetraThyraWrappers.cpp
  Thyra_get_Epetra_Operator.cpp
  )

TRILINOS_PACKAGE_ADD_LIBRARY(
  thyraepetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS thyra
  )

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
