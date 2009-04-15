
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  HbMatrix
  SOURCES cxx_main.cpp
          GenerateHbProblem.hpp
  ARGS -v
  COMM serial
  )
