

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    epetra_1DFEMNonlinear_Example
    SOURCES
    Example.C
    Basis.H
    Basis.C
    Problem_Interface.H
    Problem_Interface.C
    FiniteElementProblem.H
    FiniteElementProblem.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS -v
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )

ENDIF()
