
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    DS6_5_1
    SOURCES
    DennisSchnabel.C
    DennisSchnabel.H
    Example.C
    Problem_Interface.C
    Problem_Interface.H
    COMM serial mpi
    NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    DS6_5_1_TrustRegion
    SOURCES
    DennisSchnabel.C
    DennisSchnabel.H
    Example_TrustRegion.C
    Problem_Interface.C
    Problem_Interface.H
    COMM serial mpi
    NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

ENDIF()
