
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MultipleSolves
  SOURCES MultipleSolves.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM serial mpi
  )
