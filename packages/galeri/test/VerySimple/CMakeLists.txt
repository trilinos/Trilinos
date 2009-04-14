
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  VerySimple
  SOURCES cxx_main.cpp
  ARGS -v
  COMM serial
  )
