
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LAPACK_example
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )
