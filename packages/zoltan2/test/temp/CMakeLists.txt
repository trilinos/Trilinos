INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  paramTest
  SOURCES paramTest.cpp
  NUM_MPI_PROCS 1
  COMM mpi serial
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  multivectorTest
  SOURCES multivectorTest.cpp
  ARGS
  "401"
  NUM_MPI_PROCS 4
  COMM mpi serial
  )
