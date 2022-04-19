
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_WrappedDualView
  SOURCES
    TpetraUtils_WrappedDualView.cpp
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_WrappedDualViewMicrobench
  SOURCES
    TpetraUtils_WrappedDualViewMicrobench.cpp
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_WrappedDualViewNull
  SOURCES
    TpetraUtils_WrappedDualViewNull.cpp
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_UnitTests
  SOURCES
    TpetraUtils_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_merge
  SOURCES
    TpetraUtils_merge.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_congruent
  SOURCES
    TpetraUtils_congruent.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_gathervPrint
  SOURCES
    TpetraUtils_gathervPrint.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Behavior
  SOURCES
    behavior.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CombineMode
  SOURCES
    CombineMode.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  copyConvert
  SOURCES copyConvert.cpp
  ARGS ${ARGS}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  copyOffsets
  SOURCES copyOffsets.cpp
  ARGS ${ARGS}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  createMirrorView
  SOURCES
    createMirrorView.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  getEntryOnHost
  SOURCES
    getEntryOnHost.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test comes with its own main() function.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StaticView
  SOURCES
    StaticView.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_crsUtils
  SOURCES
    TpetraUtils_crsUtils.cpp
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  verbosePrintArray
  SOURCES
    verbosePrintArray.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE(
  createPrefix
  SOURCES
    createPrefix.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )
# In an MPI build, I'd like at least 2 processes, in order to test
# that the string prefix for verbose output includes the right rank
# and doesn't just print 0.
TRIBITS_ADD_TEST(
  createPrefix
  NAME createPrefix_MPI
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_TEST(
  createPrefix
  NAME createPrefix_no_MPI
  COMM serial
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  insertCrsIndicesThreshold
  SOURCES
    insertCrsIndicesThreshold.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
