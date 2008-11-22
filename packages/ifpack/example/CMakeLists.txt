
INCLUDE(Trilinos_Add_Executable_And_Test)

IF (${PACKAGE_NAME}_ENABLE_Galeri)
  
  TRILINOS_ADD_EXECUTABLE_AND_TEST(
    Ifpack_ex_Factory
    SOURCES Ifpack_ex_Factory.cpp
    ARGS "-v"
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )

ENDIF()
