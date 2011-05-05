
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ImportExport_UnitTests
  SOURCES
    ImportExport_UnitTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SubmapExportTests
  SOURCES
    SubmapExport
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SubmapImportTests
  SOURCES
    SubmapImport
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )
