
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Block_test_LL
  SOURCES cxx_main.cpp
  COMM serial mpi
  )

