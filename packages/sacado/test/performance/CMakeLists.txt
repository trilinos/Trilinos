INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

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

  TRIBITS_ADD_LIBRARY(
    fe_jac_fill_funcs
    HEADERS fe_jac_fill_funcs.hpp
    SOURCES fe_jac_fill_funcs.cpp
    DEPLIBS sacado
    TESTONLY
    )

  TRIBITS_ADD_EXECUTABLE(
    fad_fe_jac_fill
    SOURCES fad_fe_jac_fill.cpp
    COMM serial mpi
    XHOSTTYPE AIX
    DEPLIBS fe_jac_fill_funcs
    )

  TRIBITS_ADD_EXECUTABLE(
    fad_fe_jac_fill_range
    SOURCES fad_fe_jac_fill_range.cpp
    COMM serial mpi
    XHOSTTYPE AIX
    DEPLIBS fe_jac_fill_funcs
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

  TRIBITS_COPY_FILES_TO_BINARY_DIR(copyPythonScripts
    DEST_FILES   plot_fad_expr_depth.py)

ENDIF()
