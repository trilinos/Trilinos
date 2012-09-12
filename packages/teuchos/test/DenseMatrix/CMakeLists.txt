
INCLUDE(TribitsAddExecutableAndTest)

IF (${PACKAGE_NAME}_ENABLE_COMPLEX)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    DenseMatrix_test  
    SOURCES cxx_main.cpp
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS -v
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SymDenseMatrix_test
    SOURCES cxx_main_sym.cpp
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS -v
    STANDARD_PASS_OUTPUT
    )

ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DenseSolver_test
  SOURCES cxx_main_solver.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  ARGS -v
  STANDARD_PASS_OUTPUT
  )

