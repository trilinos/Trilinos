
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BugTests
  SOURCES
    BugTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BlankRowBugTest
  SOURCES
    BlankRowBugTest
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraBugTests
  SOURCE_FILES addA2.hb
  EXEDEPS BugTests
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CompilationTests
  SOURCES
    CompilationTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
