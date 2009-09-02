
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_MVOPTester
  SOURCES cxx_main.cpp
  COMM serial mpi
  ARGS "-v"
  STANDARD_PASS_OUTPUT
  )
