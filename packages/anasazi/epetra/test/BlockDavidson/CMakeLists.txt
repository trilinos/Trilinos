
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  BlockDavidson_solvertest
  SOURCES cxx_main_solvertest.cpp
  ARGS --verbose
  COMM serial mpi
  )
