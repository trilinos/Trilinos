
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  part_sizes
  SOURCES part_sizes.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  ARGS
)
