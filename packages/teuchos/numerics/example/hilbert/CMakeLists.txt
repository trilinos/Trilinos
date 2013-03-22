
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  hilbert
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )
