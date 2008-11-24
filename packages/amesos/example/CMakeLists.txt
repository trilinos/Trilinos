INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  example_AmesosFactory_Tridiag
  SOURCES example_AmesosFactory_Tridiag.cpp
  COMM serial mpi
  )
