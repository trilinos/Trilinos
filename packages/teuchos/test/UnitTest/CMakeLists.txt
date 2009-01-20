
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  UnitTest_UnitTests
  SOURCES
    Int_UnitTests.cpp
    TemplateFunc_UnitTests.cpp
    vector_UnitTests.cpp
    Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
PACKAGE_ADD_TEST(
  UnitTest_UnitTests
  NAME UnitTest_UnitTests_minus_one
  ARGS
    "--not-unit-test=vector_float_constAt_UnitTest"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Summary: total = 18, run = 17, passed = 17, failed = 0"
  )

PACKAGE_ADD_EXECUTABLE(
  UnitTest_BadUnitTest
  COMM serial mpi  
  SOURCES
    Int_UnitTests.cpp
    TemplateFunc_UnitTests.cpp
    Failing_UnitTest.cpp
    Teuchos_StandardUnitTestMain.cpp
  )

PACKAGE_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_final_results
  ARGS "--show-test-details=FINAL_RESULTS"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Int_BadAssignment_UnitTest"
  )
PACKAGE_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_end_result_failed
  ARGS "--teuchos-suppress-startup-banner"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "End Result: TEST FAILED"
  )
PACKAGE_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_end_result_totals
  ARGS "--teuchos-suppress-startup-banner"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Summary: total = 8, run = 8, passed = 6, failed = 2"
  )
PACKAGE_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_Int_BadAssignment_failed
  ARGS
    "--show-test-details=FINAL_RESULTS --group-name=Int --test-name=BadAssignment"
    "--show-test-details=TEST_NAMES --group-name=Int --test-name=BadAssignment"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "[FAILED]"
  )
