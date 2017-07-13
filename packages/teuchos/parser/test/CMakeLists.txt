

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Parser_UnitTests  
  SOURCES
    unit_tests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
