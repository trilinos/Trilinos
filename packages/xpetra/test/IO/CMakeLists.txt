
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IO_UnitTests
  SOURCES
    IO_UnitTests.cpp
    ../Xpetra_UnitTests.cpp
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTestsIO_cp
  SOURCE_FILES test.mtx test.mtx.bin
)
