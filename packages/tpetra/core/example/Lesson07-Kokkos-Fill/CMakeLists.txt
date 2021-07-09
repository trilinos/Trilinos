IF (Tpetra_INST_DOUBLE)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lesson07_kokkos_fill
    SOURCES 05_solve.cpp
    ARGS 
    COMM mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()
