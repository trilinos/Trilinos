
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GMRES_LL
  SOURCES cxx_main.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "All tests passed"
  )
