
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ObjectBuilder_UnitTests
  SOURCES
    ObjectBuilder_UnitTests.cpp
    Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

