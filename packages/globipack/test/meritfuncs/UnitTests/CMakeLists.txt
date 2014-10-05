INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TestLagrPolyMeritFunc1DUnitTests
  SOURCES
    TestLagrPolyMeritFunc1D_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  DEPLIBS globipack_test_meritfuncs
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
