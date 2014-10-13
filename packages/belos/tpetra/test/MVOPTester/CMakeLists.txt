

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MVOPTester_complex_test
  SOURCES cxx_main_complex.cpp
        ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )
