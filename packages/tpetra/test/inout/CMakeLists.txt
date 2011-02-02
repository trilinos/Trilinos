INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Raw_InOutTest
  SOURCES
    MatrixMarket_Raw_InOutTest.cpp
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout tpetratestlib
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_InOutTest
  SOURCES
    MatrixMarket_Tpetra_InOutTest.cpp
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout tpetratestlib
  )
