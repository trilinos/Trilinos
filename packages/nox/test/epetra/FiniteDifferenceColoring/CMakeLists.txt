

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_EpetraExt)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    FiniteDifferenceColoring
    SOURCES
    Example.C
    FiniteElementProblem.C
    FiniteElementProblem.H
    Basis.C
    Basis.H
    Problem_Interface.C
    Problem_Interface.H
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  IF (NOX_ENABLE_Isorropia)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FiniteDifferenceIsorropiaColoring
      SOURCES
      Example_IsorropiaColoring.C
      FiniteElementProblem_OwnedOnly.C
      FiniteElementProblem.H
      Basis.C
      Basis.H
      Problem_Interface.C
      Problem_Interface.H
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "Test passed!"
    )
  ENDIF()

ENDIF()
