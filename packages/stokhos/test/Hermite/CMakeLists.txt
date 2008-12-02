INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  hermite
  SOURCES hermite.cpp
  ARGS "-v"
  COMM serial mpi
  )
