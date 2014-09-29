

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA
   AND NOX_ENABLE_EpetraExt AND TPL_ENABLE_MPI)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Brusselator_XYZT_Example
    SOURCES
    Basis.H
    Basis.C
    Problem_Interface.H
    Problem_Interface.C
    Brusselator.H
    Brusselator.C
    Example.C
    COMM mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()
