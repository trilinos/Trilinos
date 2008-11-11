
INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  epetra_crsmatrix
  SOURCES cxx_main.cpp
  COMM serial mpi
  )
