
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  inout_test
  SOURCES
    cxx_main.cpp
    Poisson2dOperator.cpp
  ARGS -v
  COMM serial mpi
  )
