
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ImportExport_UnitTests
  SOURCES
    ImportExport_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Import_Union
  SOURCES
    Import_Union
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SubmapExportTests
  SOURCES
    SubmapExport
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SubmapImportTests
  SOURCES
    SubmapImport
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

IF (${PROJECT_NAME}_ENABLE_Epetra)
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

# Test whether Export to a CrsMatrix with a static graph works.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExportToStaticGraphCrsMatrix
  SOURCES ExportToStaticGraphCrsMatrix.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
  )




