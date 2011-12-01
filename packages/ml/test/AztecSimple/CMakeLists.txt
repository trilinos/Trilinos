
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AztecSimple
  SOURCES AztecSimple.c place_holder.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
