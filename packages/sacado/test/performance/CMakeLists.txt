INCLUDE(Trilinos_Add_Executable_And_Test)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/Fad)

IF (${PACKAGE_NAME}_ENABLE_Teuchos)

  TRILINOS_ADD_EXECUTABLE_AND_TEST(
    far_expr_depth
    SOURCES fad_expr_depth.cpp
    ARGS 
    COMM serial mpi
    )

  TRILINOS_ADD_EXECUTABLE_AND_TEST(
    far_expr
    SOURCES fad_expr.cpp
    ARGS 
    COMM serial mpi
    )

  TRILINOS_ADD_EXECUTABLE_AND_TEST(
    far_lj_grad
    SOURCES fad_lj_grad.cpp
    ARGS 
    COMM serial mpi
    )

  TRILINOS_ADD_EXECUTABLE_AND_TEST(
    taylor_expr
    SOURCES taylor_expr.cpp
    ARGS 
    COMM serial mpi
    )

ENDIF()