INCLUDE(Trilinos_Add_Executable_And_Test)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/Fad)

IF (${PACKAGE_NAME}_ENABLE_Teuchos)

  TRILINOS_ADD_EXECUTABLE(
    far_expr_depth
    SOURCES fad_expr_depth.cpp
    COMM serial mpi
    )

  TRILINOS_ADD_EXECUTABLE(
    far_expr
    SOURCES fad_expr.cpp
    COMM serial mpi
    )

  TRILINOS_ADD_EXECUTABLE(
    far_lj_grad
    SOURCES fad_lj_grad.cpp
    COMM serial mpi
    )

  TRILINOS_ADD_EXECUTABLE(
    taylor_expr
    SOURCES taylor_expr.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS)

    TRILINOS_ADD_TEST(
      far_expr_depth
      SOURCES fad_expr_depth.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  
    TRILINOS_ADD_TEST(
      far_expr
      SOURCES fad_expr.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  
    TRILINOS_ADD_TEST(
      far_lj_grad
      SOURCES fad_lj_grad.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  
    TRILINOS_ADD_TEST(
      taylor_expr
      SOURCES taylor_expr.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()

ENDIF()
