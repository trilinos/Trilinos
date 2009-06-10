
INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Basis.H
  FiniteElementProblem.H
  Problem_Interface.H
  Tcubed_FiniteElementProblem.H
  Pitchfork_FiniteElementProblem.H
  NormConstraint.H
  LOCALinearConstraint.H
  ) 

APPEND_SET(SOURCES
  Basis.C
  Problem_Interface.C
  Tcubed_FiniteElementProblem.C
  Pitchfork_FiniteElementProblem.C
  NormConstraint.C
  LOCALinearConstraint.C
  )

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)

  PACKAGE_ADD_LIBRARY(
    locaepetratestproblems
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    TESTONLY
    )

ENDIF()
