
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  MVOPTester
  SOURCES cxx_main.cpp 
          ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS 
  COMM serial mpi
  )
