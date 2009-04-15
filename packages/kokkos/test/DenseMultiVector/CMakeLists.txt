
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DenseMultiVector
  SOURCES cxx_main.cpp
  ARGS -v
  COMM serial
  )
