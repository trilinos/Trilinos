
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ObjectBuilder_UnitTests
  SOURCES
    ObjectBuilder_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  COMM serial mpi
  XHOSTTYPE_TEST AIX
  STANDARD_PASS_OUTPUT
  )

