INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TestLagrPolyMeritFunc1DUnitTests
  SOURCES
    TestLagrPolyMeritFunc1D_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS globipack_test_meritfuncs
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )
