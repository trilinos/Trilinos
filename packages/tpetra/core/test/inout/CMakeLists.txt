
# This test should work for any number of MPI processes, as long as
# the sparse matrix to be read in has enough rows that it may be
# distributed over all of them, and as long as Proc 0 has enough
# memory to handle all the sparse matrix's data.
IF(Tpetra_INST_INT_LONG)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_CrsMatrix_InOutTest
  SOURCES
    MatrixMarket_Tpetra_CrsMatrix_InOutTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
ENDIF()

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

# This test has an upper bound on the number of MPI processes for
# which it works, since one of the test matrices has fixed dimensions.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_Map_InOutTest
  SOURCES
    MatrixMarket_Tpetra_Map_InOutTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

# FIXME (mfh 22 Jul 2014) My recent fix to Bug 6176 breaks this test.
# Neither Andrey nor I are convinced that the test for Bug 5800 is the
# right test.  Thus, I am disabling the Bug 5800 test below, for now.

# # This test only makes sense for 1 MPI process.
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   MatrixMarket_MV_Output_Perm
#   SOURCES
#     MatrixMarket_MV_Output_Perm.cpp
#     ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   ARGS
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
#   )

# This test doesn't actually run as part of the check-in tests; it
# just builds.  It's good for testing Matrix Market input and output
# for Map and CrsMatrix, with arbitrary files.  The test should work
# with or without MPI support.
IF (Tpetra_INST_INT_LONG)
TRIBITS_ADD_EXECUTABLE(
  MatrixMarket_Tpetra_CrsMatrix_FileTest
  SOURCES
    MatrixMarket_Tpetra_CrsMatrix_FileTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )
ENDIF()

