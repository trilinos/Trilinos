
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ImportExport_UnitTests
  SOURCES
    ImportExport_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense in an MPI build with > 1 processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ImportExport_ImportConstructExpert
  SOURCES
    ImportExport_ImportConstructExpert.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2-4
  STANDARD_PASS_OUTPUT
  )



### Import memory testing for GDSW ###
# This test only makes sense in an MPI build with > 1 processes.
TRIBITS_ADD_EXECUTABLE(
  ImportExport_GDSWStyle
  SOURCES
  ImportExport_GDSWStyle.cpp
  COMM mpi
  CATEGORIES BASIC PERFORMANCE
  )

TRIBITS_ADD_TEST(
  ImportExport_GDSWStyle
  NAME GDSWStyle_BaselineTpetra
  ARGS "--case=1 --watchr-output"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  RUN_SERIAL
  CATEGORIES BASIC PERFORMANCE
  )

TRIBITS_ADD_TEST(
  ImportExport_GDSWStyle
  NAME GDSWStyle_ImportSquareMatrix
  ARGS "--case=2 --watchr-output"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  RUN_SERIAL
  CATEGORIES BASIC PERFORMANCE
  )



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Import_Union
  SOURCES
    Import_Union.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  ARGS "--globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug9116
  SOURCES
    Bug9116.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SubmapExportTests
  SOURCES
    SubmapExport.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SubmapImportTests
  SOURCES
    SubmapImport.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnpackLongRows
    SOURCES
      UnpackLongRows.cpp
    COMM mpi
    NUM_MPI_PROCS 4
    STANDARD_PASS_OUTPUT
  )

IF (Tpetra_ENABLE_DEPRECATED_CODE AND ${PROJECT_NAME}_ENABLE_Epetra)
  IF(NOT Trilinos_NO_32BIT_GLOBAL_INDICES AND Tpetra_INST_INT_INT)
    # Tpetra bug 5430:
    # https://software.sandia.gov/bugzilla/show_bug.cgi?id=5430
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
ENDIF ()

IF (Tpetra_INST_INT_LONG_LONG OR Tpetra_INST_INT_LONG OR Tpetra_INST_INT_INT OR Tpetra_INST_INT_UNSIGNED_LONG OR Tpetra_INST_INT_UNSIGNED)
  # Test whether Export to a CrsMatrix with a static graph works.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ExportToStaticGraphCrsMatrix
    SOURCES
      ExportToStaticGraphCrsMatrix.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 4
    )
ENDIF()

# Test Tpetra::Details::makeOptimizedColMap and
# Tpetra::Details::makeOptimizedColMapAndImport.
#
# This test really only makes sense with an MPI build, since
# optimizing the column Map involves reordering remote indices.  The
# test already exercises the case where there are no remote indices.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MakeOptColMap
  SOURCES
    MakeOptColMap.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
  )

# mfh 18 Jul 2014: Bug 6170 test.  Thanks to Clark Dohrmann for
# contributing this test.  I modified it slightly to avoid file I/O
# and to catch and print exception messages on all processes.  The
# test must be run with exactly 8 MPI processes.
#
# NOTE: The default behavior of the Trilinos test framework is to
# exclude tests that require more than 4 MPI processes.  Thus, this
# test will NOT run by default.  To change this, set the CMake option
# MPI_EXEC_MAX_NUMPROCS to 8 or more in your Trilinos configuration.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug6170
  SOURCES
    Bug6170.cpp
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 8
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ReverseCommunication_Issue227
  SOURCES
    ReverseCommunication_Issue227.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

# mfh 17 May 2016: Test for Github Issue #114.
# Test originally written by Chris Siefert.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue_114
  SOURCES Issue_114.cpp
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
  )

# mfh 13 Sep 2016: Test for Github Issue #607.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue_607
  SOURCES
    Issue_607.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1-4
  )

# mfh 01 Feb 2018: Test for Github Issue #2198.
# It should always build, but it's only useful
# to run in an MPI build.
TRIBITS_ADD_EXECUTABLE(
  Issue_2198
  SOURCES
    Issue_2198.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

# mfh 10 May 2018: See GitHub Issue #2564 for why we split the GitHub
# Issue #2198 test into two separate cases: the first needing 3 MPI
# processes, and the second needing 5 MPI processes.

TRIBITS_ADD_TEST(
  Issue_2198
  NAME Issue_2198_no_more_than_4_processes
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Issue_2198
  NAME Issue_2198_exactly_5_processes
  COMM mpi
  NUM_MPI_PROCS 5
  STANDARD_PASS_OUTPUT
  )

# Always build the test, whether or not MPI is enabled,
# but only run it when MPI is enabled, because it is only
# meaningful in that case.
TRIBITS_ADD_EXECUTABLE(
  Issue_3968
  SOURCES
    Issue_3968.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  Issue_3968
  NAME Issue_3968_mpi
  ARGS "--globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )
