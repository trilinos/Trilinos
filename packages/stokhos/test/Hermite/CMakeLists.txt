INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  hermite
  SOURCES hermite.cpp
  ARGS "-v"
  COMM serial mpi
  )
