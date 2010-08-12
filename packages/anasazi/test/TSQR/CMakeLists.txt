INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   TSQR_Verify
#   SOURCES Tsqr_Tester.cpp
#   COMM serial mpi
#   ARGS 
#   STANDARD_PASS_OUTPUT
# )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_CombineVerify
  SOURCES Tsqr_TestCombine.cpp
  COMM serial
  ARGS "--verify --nrows=1000 --ncols=5"
  STANDARD_PASS_OUTPUT
)

