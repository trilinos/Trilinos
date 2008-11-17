
INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  AztecOO_test
  SOURCES cxx_main.cpp
  ARGS -v
  COMM serial mpi
  )
