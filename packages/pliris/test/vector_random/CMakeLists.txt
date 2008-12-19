
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  vector_random
  SOURCES cxx_main.cpp
  ARGS 1000 1
  COMM serial mpi
  )
