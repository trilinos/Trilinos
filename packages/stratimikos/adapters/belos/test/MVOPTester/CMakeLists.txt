INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MVOPTesterThyra_test
  SOURCES cxx_main.cpp
  ARGS "-v"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ThyraTsqrAdaptor_test
  SOURCES ThyraTsqrAdaptorTester.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
