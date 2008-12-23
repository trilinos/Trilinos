
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Permutation_test
  SOURCES cxx_main.cpp
  COMM serial mpi
  )
