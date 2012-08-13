INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  paramTest
  SOURCES paramTest.cpp
  NUM_MPI_PROCS 4
  COMM mpi serial
  )
