INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTestHarness_Parallel_UnitTests
  SOURCES UnitTestHarness_Parallel_UnitTests.cpp
    ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_ADVANCED_TEST(
  UnitTestHarness_Parallel_UnitTests_MPI
    TEST_0
      EXEC UnitTestHarness_Parallel_UnitTests
      PASS_REGULAR_EXPRESSION_ALL
        "NOTE: Global reduction shows failures on other processes"
        ".*FAILED.* UnitTestHarness_nonRankFails_UnitTest"
        ".*FAILED.* UnitTestHarness_nonRankThrowsTeuchosExcept_UnitTest"
        ".*FAILED.* UnitTestHarness_nonRankThrowsIntExcept_UnitTest"
  COMM mpi
  )
