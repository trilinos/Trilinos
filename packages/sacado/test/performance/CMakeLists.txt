INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/Fad)
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test/utils)

IF (${PACKAGE_NAME}_ENABLE_Teuchos)

  TRIBITS_ADD_EXECUTABLE(
    fad_expr_depth
    SOURCES fad_expr_depth.cpp fad_expr_funcs.hpp fad_expr_funcs.cpp
    COMM serial mpi
    XHOSTTYPE AIX
    )

  TRIBITS_ADD_EXECUTABLE(
    fad_expr
    SOURCES fad_expr.cpp
    COMM serial mpi
    XHOSTTYPE AIX
    )

  TRIBITS_ADD_EXECUTABLE(
    fad_lj_grad
    SOURCES fad_lj_grad.cpp
    COMM serial mpi
    XHOSTTYPE AIX
    )

  TRIBITS_ADD_EXECUTABLE(
    taylor_expr
    SOURCES taylor_expr.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    fad_fe_jac_fill
    SOURCES fad_fe_jac_fill.cpp
    COMM serial mpi
    XHOSTTYPE AIX
    )

  TRIBITS_ADD_EXECUTABLE(
    rad_fe_adj_fill
    SOURCES rad_fe_adj_fill.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    fad_blas
    SOURCES fad_blas.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS)

    TRIBITS_ADD_TEST(
      fad_expr_depth
      SOURCES fad_expr_depth.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      XHOSTTYPE AIX
      )
  
    TRIBITS_ADD_TEST(
      fad_expr
      SOURCES fad_expr.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      XHOSTTYPE AIX
      )
  
    TRIBITS_ADD_TEST(
      fad_lj_grad
      SOURCES fad_lj_grad.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      XHOSTTYPE AIX
      )
  
    TRIBITS_ADD_TEST(
      taylor_expr
      SOURCES taylor_expr.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_TEST(
      fad_fe_jac_fill
      SOURCES fad_fe_jac_fill.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      XHOSTTYPE AIX
      )

    TRIBITS_ADD_TEST(
      rad_fe_adj_fill
      SOURCES rad_fe_adj_fill.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_TEST(
      fad_blas
      SOURCES fad_blas.cpp
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()

ENDIF()
