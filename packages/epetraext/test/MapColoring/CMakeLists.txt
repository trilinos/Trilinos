
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EpetraExt_MapColoring_test
  SOURCES cxx_main.cpp
  COMM serial
  )
