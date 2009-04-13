INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_probing
  SOURCES probing.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  ARGS
  )
