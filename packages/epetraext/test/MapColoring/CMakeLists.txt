
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  EpetraExt_MapColoring_test
  SOURCES cxx_main.cpp
  COMM serial
  )
