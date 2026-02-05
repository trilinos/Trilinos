
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MultiVector_Adapter_UnitTests
  SOURCES
    Tpetra_MultiVector_Adapter_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_CrsMatrix_Adapter_UnitTests
  SOURCES
    Tpetra_CrsMatrix_Adapter_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

