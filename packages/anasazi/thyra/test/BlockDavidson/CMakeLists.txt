
INCLUDE(Trilinos_Add_Executable_And_Test)

IF(Anasazi_ENABLE_Epetra)
  TRILINOS_ADD_EXECUTABLE_AND_TEST(
    BlockDavidsonThyra_test
    SOURCES cxx_main.cpp
    ARGS 
    COMM serial mpi
    )
ENDIF()
