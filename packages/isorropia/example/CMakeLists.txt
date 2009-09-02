
INCLUDE(PackageAddExecutableAndTest)


ADD_SUBDIRECTORY(geometric)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  matrix_1
  SOURCES matrix_1.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
