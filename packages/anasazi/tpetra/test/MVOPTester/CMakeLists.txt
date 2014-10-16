

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MVOPTester
  SOURCES cxx_main.cpp 
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVecTraitsTest2
  SOURCES
    MultiVecTraitsTest2.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )
