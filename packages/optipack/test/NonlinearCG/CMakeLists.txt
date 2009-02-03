INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  NonlinearCGUnitTests
  SOURCES
    NonlinearCG_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS optipack_test_models
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )
