INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

# This test should work for any number of MPI processes, as long as
# the sparse matrix to be read in has enough rows that it may be
# distributed over all of them.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_InOutTest
  SOURCES
    MatrixMarket_Tpetra_InOutTest.cpp
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test should work for any number of MPI processes, as long as
# the dense matrix to be read in has enough rows that it may be
# distributed over all of them.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_MultiVector_InOutTest
  SOURCES
    MatrixMarket_Tpetra_MultiVector_InOutTest.cpp
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

# This test should work for any number of MPI processes, as long as
# Proc 0 has enough memory to handle all the Map data.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_Map_InOutTest
  SOURCES
    MatrixMarket_Tpetra_Map_InOutTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
  
