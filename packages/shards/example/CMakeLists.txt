
INCLUDE(Trilinos_Add_Executable_And_Test)

IF(Shards_ENABLE_Teuchos)

  TRILINOS_ADD_EXECUTABLE_AND_TEST(
    example_01
    SOURCES example_01.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

ENDIF()