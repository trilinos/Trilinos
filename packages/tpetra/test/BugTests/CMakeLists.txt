
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BugTests
  SOURCES
    BugTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlankRowBugTest
  SOURCES
    BlankRowBugTest
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

# Tpetra bug 5401:
# https://software.sandia.gov/bugzilla/show_bug.cgi?id=5401
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NegativeBaseIndexTest
  SOURCES
    NegativeBaseIndexTest
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraBugTests
  SOURCE_FILES addA2.hb
  EXEDEPS BugTests
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CompilationTests
  SOURCES
    CompilationTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

IF (${PROJECT_NAME}_ENABLE_Epetra)
  # Tpetra bug 5430:
  # https://software.sandia.gov/bugzilla/show_bug.cgi?id=5401
  #
  # The bug only manifests with 4 MPI processes.  The test currently
  # requires Epetra; it compares Epetra and Tpetra results.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ImportBug5430
    SOURCES
    ImportBug5430Test.cpp
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 4
    )
ENDIF ()
