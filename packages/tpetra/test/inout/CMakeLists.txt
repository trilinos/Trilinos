INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# FIXME (mfh 02 Feb 2011) I don't think this test works with > 1 MPI
# process.  I should fix it so it does.  For now, I just require that
# the test only use 1 MPI process.
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Raw_InOutTest
  SOURCES
    MatrixMarket_Raw_InOutTest.cpp
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout tpetratestlib
  )

# This test should work for any number of MPI processes, as long as
# the sparse matrix to be read in has enough rows that it may be
# distributed over all of them.
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_InOutTest
  SOURCES
    MatrixMarket_Tpetra_InOutTest.cpp
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetra tpetrainout tpetratestlib
  )
