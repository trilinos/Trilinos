
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Test_KLU
  SOURCES cxx_main.cpp
  ARGS 
  COMM serial
  )
