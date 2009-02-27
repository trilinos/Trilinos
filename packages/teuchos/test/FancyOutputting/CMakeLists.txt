
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  FancyOutputting_test  
  SOURCES FancyOutputting_test.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
