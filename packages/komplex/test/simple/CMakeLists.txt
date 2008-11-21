
INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  simple
  SOURCES simple.c
  ARGS -v
  COMM serial mpi
  )
