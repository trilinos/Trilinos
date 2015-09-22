

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_PETSC)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    petsc_ds6.5.1_example
    SOURCES Example.C DennisSchnabel.H DennisSchnabel.C Problem_Interface.H Problem_Interface.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS "-nox_linesearch_type polynomial"
    PASS_REGULAR_EXPRESSION "Converged!"
    )

ENDIF()
