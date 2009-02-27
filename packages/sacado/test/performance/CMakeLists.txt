INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/Fad)

IF (${PACKAGE_NAME}_ENABLE_Teuchos)

  PACKAGE_ADD_EXECUTABLE(
    fad_expr_depth
    SOURCES fad_expr_depth.cpp
    DEPLIBS sacadoutils
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE(
    fad_expr
    SOURCES fad_expr.cpp
    DEPLIBS sacadoutils
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE(
    fad_lj_grad
    SOURCES fad_lj_grad.cpp
    DEPLIBS sacadoutils
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE(
    taylor_expr
    SOURCES taylor_expr.cpp
    DEPLIBS sacadoutils
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE(
    fad_fe_fill
    SOURCES fad_fe_fill.cpp
    DEPLIBS sacadoutils
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS)

    PACKAGE_ADD_TEST(
      fad_expr_depth
      SOURCES fad_expr_depth.cpp
      DEPLIBS sacadoutils
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  
    PACKAGE_ADD_TEST(
      fad_expr
      SOURCES fad_expr.cpp
      DEPLIBS sacadoutils
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  
    PACKAGE_ADD_TEST(
      fad_lj_grad
      SOURCES fad_lj_grad.cpp
      DEPLIBS sacadoutils
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  
    PACKAGE_ADD_TEST(
      taylor_expr
      SOURCES taylor_expr.cpp
      DEPLIBS sacadoutils
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

    PACKAGE_ADD_TEST(
      fad_fe_fill
      SOURCES fad_fe_fill.cpp
      DEPLIBS sacadoutils
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()

ENDIF()
