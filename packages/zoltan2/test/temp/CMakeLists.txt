INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  testTrilinos 
  SOURCES testTrilinos.cpp
  NUM_MPI_PROCS 4
  COMM mpi serial
  )
