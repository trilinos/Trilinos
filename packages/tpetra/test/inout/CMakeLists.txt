INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

# This test should work for any number of MPI processes, as long as
# the sparse matrix to be read in has enough rows that it may be
# distributed over all of them, and as long as Proc 0 has enough
# memory to handle all the sparse matrix's data.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_CrsMatrix_InOutTest
  SOURCES
    MatrixMarket_Tpetra_CrsMatrix_InOutTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test should work for any number of MPI processes, as long as
# the dense matrix to be read in has enough rows that it may be
# distributed over all of them, and as long as Proc 0 has enough
# memory to handle all the dense matrix's data.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_MultiVector_InOutTest
  SOURCES
    MatrixMarket_Tpetra_MultiVector_InOutTest.cpp
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test should work for any number of MPI processes, as long as
# Proc 0 has enough memory to handle all the Map's data.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_Map_InOutTest
  SOURCES
    MatrixMarket_Tpetra_Map_InOutTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense for 1 MPI process.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_MV_Output_Perm
  SOURCES
    MatrixMarket_MV_Output_Perm.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test doesn't actually run as part of the check-in tests; it
# just builds.  It's good for testing Matrix Market input and output
# for Map and CrsMatrix, with arbitrary files.  The test should work
# with or without MPI support.
TRIBITS_ADD_EXECUTABLE(
  MatrixMarket_Tpetra_CrsMatrix_FileTest
  SOURCES
    MatrixMarket_Tpetra_CrsMatrix_FileTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )
  
