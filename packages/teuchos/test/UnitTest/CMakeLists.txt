
INCLUDE(TribitsAddExecutableAndTest)


#
# UnitTest_UnitTests
#


TRIBITS_ADD_EXECUTABLE(
  UnitTest_UnitTests
  SOURCES
    UnitTest_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )


TRIBITS_ADD_TEST(
  UnitTest_UnitTests
  POSTFIX_AND_ARGS_0 test_names --details=TEST_NAMES
  POSTFIX_AND_ARGS_1 final_results --details=FINAL_RESULTS
  POSTFIX_AND_ARGS_2 all --details=ALL
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


#
# Simple_UnitTests
#

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Simple_UnitTests
  SOURCES
    Int_UnitTests.cpp
    TemplateFunc_UnitTests.cpp
    vector_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_COMPLEX AND ${PACKAGE_NAME}_ENABLE_FLOAT)
  SET(PASS_STRING "Summary: total = 18, run = 17, passed = 17, failed = 0")
ELSEIF (${PACKAGE_NAME}_ENABLE_COMPLEX AND NOT ${PACKAGE_NAME}_ENABLE_FLOAT)
  SET(PASS_STRING "Summary: total = 16, run = 15, passed = 15, failed = 0")
ELSEIF (NOT ${PACKAGE_NAME}_ENABLE_COMPLEX AND ${PACKAGE_NAME}_ENABLE_FLOAT)
  SET(PASS_STRING "Summary: total = 16, run = 15, passed = 15, failed = 0")
ELSEIF (NOT ${PACKAGE_NAME}_ENABLE_COMPLEX AND NOT ${PACKAGE_NAME}_ENABLE_FLOAT)
  SET(PASS_STRING "Summary: total = 15, run = 14, passed = 14, failed = 0")
ENDIF()

TRIBITS_ADD_TEST(
  Simple_UnitTests
  NAME_POSTFIX minus_one
  ARGS "--not-unit-test=vector_float_constAt_UnitTest"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "${PASS_STRING}"
  )


#
# UnitTest_TimeMonitor_UnitTests
#

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


#
# UnitTest_BadUnitTest
#

TRIBITS_ADD_EXECUTABLE(
  UnitTest_BadUnitTest
  COMM serial mpi  
  SOURCES
    Int_UnitTests.cpp
    TemplateFunc_UnitTests.cpp
    Failing_UnitTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )

TRIBITS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_final_results
  ARGS "--show-test-details=FINAL_RESULTS"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Int_BadAssignment_UnitTest"
  )

TRIBITS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_end_result_failed
  ARGS "--teuchos-suppress-startup-banner"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "End Result: TEST FAILED"
  )

TRIBITS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_end_result_totals
  ARGS "--teuchos-suppress-startup-banner"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION ", failed = 2"
  )

TRIBITS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_Int_BadAssignment_failed
  ARGS
    "--show-test-details=FINAL_RESULTS --group-name=Int --test-name=BadAssignment"
    "--show-test-details=TEST_NAMES --group-name=Int --test-name=BadAssignment"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "[FAILED]"
  )
