

TRIBITS_ADD_EXECUTABLE(
  UnitTest_TimeMonitor_UnitTests
  SOURCES
    TimeMonitor_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )

TRIBITS_ADD_ADVANCED_TEST(
  UnitTest_TimeMonitor_UnitTests
  TEST_0 EXEC UnitTest_TimeMonitor_UnitTests
    ARGS --test=someFunction_timed --show-timer-summary
    PASS_REGULAR_EXPRESSION_ALL
      "someFunction1"
      "someFunction2"
  OVERALL_NUM_MPI_PROCS 1
  XHOSTTYPE AIX
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue3868
  SOURCES
    TimeMonitor_Issue3868.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SyncTimeMonitor
  SOURCES
    SyncTimeMonitor.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 4
  )
