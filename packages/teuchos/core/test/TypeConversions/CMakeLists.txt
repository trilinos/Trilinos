
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TypeConversions_UnitTest
  SOURCES TypeConversions_UnitTest.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
