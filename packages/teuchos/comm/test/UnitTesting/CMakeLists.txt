INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)


INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTestHarness_Parallel_UnitTests
  SOURCES UnitTestHarness_Parallel_UnitTests.cpp
    ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  ARGS "--details=ALL"
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_ADVANCED_TEST(
  UnitTestHarness_Parallel_UnitTests_MPI_2
  TEST_0
    EXEC UnitTestHarness_Parallel_UnitTests
    MESSAGE "Show default output on proc 0"
    ARGS --output-show-proc-rank
    PASS_REGULAR_EXPRESSION_ALL
      "NOTE: Global reduction shows failures on other processes"
      "NOTE: Unit test failed on processes = {1}"
      ".*FAILED.* UnitTestHarness_nonRootFails_UnitTest"
      ".*FAILED.* UnitTestHarness_nonRootThrowsTeuchosExcept_UnitTest"
      ".*FAILED.* UnitTestHarness_nonRootThrowsIntExcept_UnitTest"
    NUM_MPI_PROCS 2
  TEST_1
    EXEC UnitTestHarness_Parallel_UnitTests
    MESSAGE "Show output on proc 1"
    ARGS --output-show-proc-rank --output-to-root-rank-only=1
    PASS_REGULAR_EXPRESSION_ALL
      "NOTE: Unit test failed on processes = {1}"
      ".*FAILED.* UnitTestHarness_nonRootFails_UnitTest"
      ".*FAILED.* UnitTestHarness_nonRootThrowsTeuchosExcept_UnitTest"
      ".*FAILED.* UnitTestHarness_nonRootThrowsIntExcept_UnitTest"
    NUM_MPI_PROCS 2
  COMM mpi
  )

TRIBITS_ADD_ADVANCED_TEST(
  UnitTestHarness_Parallel_UnitTests_MPI_4
  TEST_0
    EXEC UnitTestHarness_Parallel_UnitTests
    MESSAGE "Show default output on proc 0"
    ARGS --output-show-proc-rank
    PASS_REGULAR_EXPRESSION_ALL
      "NOTE: Global reduction shows failures on other processes"
      "NOTE: Unit test failed on processes = {1, 3}"
      ".*FAILED.* UnitTestHarness_nonRootFails_UnitTest"
      ".*FAILED.* UnitTestHarness_nonRootThrowsTeuchosExcept_UnitTest"
      ".*FAILED.* UnitTestHarness_nonRootThrowsIntExcept_UnitTest"
    NUM_MPI_PROCS 4
  TEST_1
    EXEC UnitTestHarness_Parallel_UnitTests
    MESSAGE "Show output on proc 1"
    ARGS --output-show-proc-rank --output-to-root-rank-only=1
    PASS_REGULAR_EXPRESSION_ALL
      "NOTE: Unit test failed on processes = {1, 3}"
      ".*FAILED.* UnitTestHarness_nonRootFails_UnitTest"
      ".*FAILED.* UnitTestHarness_nonRootThrowsTeuchosExcept_UnitTest"
      ".*FAILED.* UnitTestHarness_nonRootThrowsIntExcept_UnitTest"
    NUM_MPI_PROCS 4
  COMM mpi
  )

TRIBITS_ADD_ADVANCED_TEST(
  UnitTestHarness_Parallel_UnitTests_MPI_no_reduce
    TEST_0
      EXEC UnitTestHarness_Parallel_UnitTests
      ARGS --test=nonRootFails --no-globally-reduce-test-result
      PASS_REGULAR_EXPRESSION_ALL
        "UnitTestHarness_nonRootFails_UnitTest ... .Passed."
      NUM_MPI_PROCS 2
  COMM mpi
  )
