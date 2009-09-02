
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TwoSquares
  SOURCES TwoSquares.cpp
  ARGS 
  COMM serial
  )
