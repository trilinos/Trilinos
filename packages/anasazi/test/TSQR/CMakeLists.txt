INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_Verify
  SOURCES Tsqr_Tester.cpp
  COMM serial mpi
  ARGS 
    "--verbose"
    "--verify"
#    "--debug"
  STANDARD_PASS_OUTPUT
)
