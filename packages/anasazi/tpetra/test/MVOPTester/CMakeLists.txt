
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MVOPTester
  SOURCES cxx_main.cpp 
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )
