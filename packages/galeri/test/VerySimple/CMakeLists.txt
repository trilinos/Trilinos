
INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  VerySimple
  SOURCES cxx_main.cpp
  ARGS 
  COMM serial mpi
  )
