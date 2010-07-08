INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MultiVector_Adapter_UnitTests
  SOURCES
    Tpetra_MultiVector_Adapter_UnitTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Tpetra_CrsMatrix_Adapter_UnitTests
  SOURCES
    Tpetra_CrsMatrix_Adapter_UnitTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
