INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_example
  SOURCES simple.c
  ARGS -v
  COMM serial mpi
  )
