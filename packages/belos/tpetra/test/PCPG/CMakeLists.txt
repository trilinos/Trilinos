TRIBITS_ADD_EXECUTABLE_AND_TEST(
  pcpg_tpetraex
  SOURCES test_pcpg_tpetraex.cpp
  COMM mpi
  NUM_MPI_PROCS 4
  )
