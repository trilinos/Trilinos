
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  inout_test
  SOURCES
    cxx_main.cpp
    Poisson2dOperator.cpp
  COMM serial mpi
  )
