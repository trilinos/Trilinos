TRIBITS_ADD_EXECUTABLE(
  ImportExport2_UnitTests
  SOURCES
    ImportExport2_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  ImportExport2_UnitTests
  NAME ImportExport2_UnitTests_Send
  COMM serial mpi
  ARGS "--distributor-send-type=Send --globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  ImportExport2_UnitTests
  NAME ImportExport2_UnitTests_ISend
  COMM serial mpi
  ARGS "--distributor-send-type=Isend --globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AsyncTransfer_UnitTests
  SOURCES
    AsyncTransfer_UnitTests.cpp
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
