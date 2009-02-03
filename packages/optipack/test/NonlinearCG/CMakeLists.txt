INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE(
  NonlinearCGUnitTests
  SOURCES
    NonlinearCG_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS optipack_test_models
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_TEST(
  NonlinearCGUnitTests
  NAME NonlinearCGUnitTests_serial
  ARGS
    " "
    "--local-dim=5"
    "--local-dim=10"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_TEST(
  NonlinearCGUnitTests
  NAME NonlinearCGUnitTests_mpi
  ARGS
    "--local-dim=4"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )
