
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MVOPTester_complex_test
  SOURCES cxx_main_complex.cpp
        ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS 
  COMM serial mpi
  )
