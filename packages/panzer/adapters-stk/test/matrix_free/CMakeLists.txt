TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER ${PANZER_UNIT_TEST_MAIN})

# Belos CG preconditioned by MueLu p-multigrid with Chebyshev smoothing (stratimikos-chebyshev.xml), and by Ifpack2
# Jacobi (stratimikos-jacobi.xml), on the assembled matrix and on the matrix-free operator
TRIBITS_ADD_EXECUTABLE(
  Test_Matrix_Free
  SOURCES test_fe_matrix_free.cpp
  COMM mpi
)

TRIBITS_ADD_TEST(
  Test_Matrix_Free
  NAME_POSTFIX MueLuChebyshev
  ARGS "--stratimikosParams=stratimikos-chebyshev.xml"
  NUM_MPI_PROCS 2
  COMM mpi
  PASS_REGULAR_EXPRESSION "TEST PASSED"
)

TRIBITS_ADD_TEST(
  Test_Matrix_Free
  NAME_POSTFIX Jacobi
  ARGS "--stratimikosParams=stratimikos-jacobi.xml"
  NUM_MPI_PROCS 2
  COMM mpi
  PASS_REGULAR_EXPRESSION "TEST PASSED"
)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/stratimikos-chebyshev.xml"
               "${CMAKE_CURRENT_BINARY_DIR}/stratimikos-chebyshev.xml" COPYONLY)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/stratimikos-jacobi.xml"
               "${CMAKE_CURRENT_BINARY_DIR}/stratimikos-jacobi.xml" COPYONLY)
