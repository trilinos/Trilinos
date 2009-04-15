
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  OskiMatrix
  SOURCES cxx_main.cpp
          GenerateOskiProblem.hpp
  ARGS -v
  COMM serial
  )
