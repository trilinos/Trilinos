
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Transpose_test
  SOURCES cxx_main.cpp
  COMM serial mpi
  )
