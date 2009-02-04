
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  StatusTest_test
  SOURCES cxx_main.cpp
  COMM serial mpi
  ARGS
    "--verbose"
#    "--debug"
  STANDARD_PASS_OUTPUT
  )
