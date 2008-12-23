
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BaseSparseMultiply
  SOURCES cxx_main.cpp
          scsrmm2.f
          GenerateHbProblem.hpp
  ARGS -v
  COMM serial
  )
