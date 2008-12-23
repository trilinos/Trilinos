INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  example_AmesosFactory_Tridiag
  SOURCES example_AmesosFactory_Tridiag.cpp
  COMM serial mpi
  )
