INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  test_epetra_adapters
  SOURCES test_epetra_adapters.cpp
  ARGS --verbose
  COMM serial mpi
  )
