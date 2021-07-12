TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ImportExport2_UnitTests
  SOURCES
    ImportExport2_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  ARGS "--globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ImportExport2_CrsSortingUtils
  SOURCES
    ImportExport2_CrsSortingUtils.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
