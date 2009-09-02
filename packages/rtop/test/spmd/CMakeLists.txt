INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  runSpmdTests
  SOURCES runSpmdTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
