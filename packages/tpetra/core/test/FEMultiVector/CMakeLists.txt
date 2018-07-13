TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEMultiVector_UnitTests
  SOURCES
    FEMultiVector_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  ARGS "--globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEMultiVector_Fix3101
  SOURCES
    Fix3101.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

