INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  RTIExample
  SOURCES RTIExample
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)
