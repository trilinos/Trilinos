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

TRIBITS_ADD_TEST(
  ImportExport2_UnitTests
  NAME ImportExport2_UnitTests_Alltoall
  COMM serial mpi
  ARGS "--distributor-send-type=Alltoall --globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  STANDARD_PASS_OUTPUT
  )

if (TPL_ENABLE_MPI)
TRIBITS_ADD_TEST(
  ImportExport2_UnitTests
  NAME ImportExport2_UnitTests_Ialltofewv
  COMM mpi
  ARGS "--distributor-send-type=Ialltofewv --globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  STANDARD_PASS_OUTPUT
)
endif()

IF (${PACKAGE_NAME}_ENABLE_mpi_advance)
  message(STATUS "Tpetra: adding MPI Advance ImportExport2 unit tests because ${PACKAGE_NAME}_ENABLE_mpi_advance is truthy")

  TRIBITS_ADD_TEST(
    ImportExport2_UnitTests
    NAME ImportExport2_UnitTests_MpiAdvanceAlltoall
    COMM serial mpi
    ARGS "--distributor-send-type=MpiAdvanceAlltoall --globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
    STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_TEST(
      ImportExport2_UnitTests
      NAME ImportExport2_UnitTests_MpiAdvanceNbralltoallv
      COMM serial mpi
      ARGS "--distributor-send-type=MpiAdvanceNbralltoallv --globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
      STANDARD_PASS_OUTPUT
      )
endif()

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
