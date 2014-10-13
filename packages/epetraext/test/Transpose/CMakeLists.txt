

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Transpose_test
  SOURCES cxx_main.cpp
  ARGS -v
  COMM serial mpi
  )
