INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES hermite_example.cpp
  ARGS -v
  COMM serial mpi
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  product_basis_example
  SOURCES product_basis_example.cpp
  ARGS -v
  COMM serial mpi
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  dakota_pce
  SOURCES dakota_pce.cpp
  ARGS -v
  COMM serial mpi
  )
