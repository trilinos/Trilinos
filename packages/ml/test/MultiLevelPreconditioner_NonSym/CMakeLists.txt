
INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  MultiLevelPreconditioner_NonSym
  SOURCES MultiLevelPreconditioner_NonSym.cpp
  ARGS 
  COMM serial
  )
