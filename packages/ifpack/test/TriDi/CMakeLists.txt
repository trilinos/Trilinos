

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TriDi
  SOURCES tridi_main.cpp
  ARGS "-v"
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

