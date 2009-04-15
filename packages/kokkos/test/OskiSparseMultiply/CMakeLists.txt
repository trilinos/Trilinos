
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  OskiSparseMultiply
  SOURCES cxx_main.cpp
          GenerateOskiProblem.hpp
          ../BaseSparseMultiply/scsrmm2.f
  ARGS -v
  COMM serial
  )
