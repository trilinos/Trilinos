
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  iallreduce
  SOURCES
    iallreduce
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
