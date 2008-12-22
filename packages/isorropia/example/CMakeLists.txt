
INCLUDE(Trilinos_Add_Executable_And_Test)


ADD_SUBDIRECTORY(geometric)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  matrix_1
  SOURCES matrix_1.cpp
  DEPLIBS isorropia_test_utils
  COMM serial mpi
  )
