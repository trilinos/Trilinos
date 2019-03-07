TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/Fad)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test/utils)

# Don't build these tests when Cuda is enabled -- they take forever to compile
# using the nvcc wrapper
IF (${PACKAGE_NAME}_ENABLE_TeuchosCore AND NOT Kokkos_ENABLE_Cuda)

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
    TESTONLY
    )

  TRIBITS_ADD_EXECUTABLE(
    fad_fe_jac_fill
    SOURCES fad_fe_jac_fill.cpp
    COMM serial mpi
    XHOSTTYPE AIX
    TESTONLYLIBS fe_jac_fill_funcs
    )

  TRIBITS_ADD_EXECUTABLE(
    fad_fe_jac_fill_range
    SOURCES fad_fe_jac_fill_range.cpp
    COMM serial mpi
    XHOSTTYPE AIX
    TESTONLYLIBS fe_jac_fill_funcs
    )

  TRIBITS_ADD_EXECUTABLE(
    rad_fe_adj_fill
    SOURCES rad_fe_adj_fill.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_TeuchosNumerics)
    TRIBITS_ADD_EXECUTABLE(
      fad_blas
      SOURCES fad_blas.cpp
      COMM serial mpi
      )
  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(copyPythonScripts
    DEST_FILES   plot_fad_expr_depth.py)

ENDIF()

IF (Sacado_ENABLE_TeuchosCore AND Sacado_ENABLE_KokkosCore)

  TRIBITS_ADD_EXECUTABLE(
    fad_kokkos_view
    SOURCES fad_kokkos_view.cpp
    COMM serial mpi
    )

  # These tests do not compile with gcc 4.7.x because it doesn't properly
  # support lambdas.  See github issue #854
  IF(NOT ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")))
    TRIBITS_ADD_EXECUTABLE(
      fad_kokkos_hierarchical
      SOURCES fad_kokkos_hierarchical.cpp
      COMM serial mpi
      )

#    TRIBITS_ADD_EXECUTABLE(
#      fad_kokkos_mat_vec_perf
#      SOURCES fad_kokkos_mat_vec_perf.cpp
#      COMM serial mpi
#      )
  ENDIF()

ENDIF()

ADD_SUBDIRECTORY(fenl_assembly)
ADD_SUBDIRECTORY(fenl_assembly_view)
