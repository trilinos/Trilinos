TRIBITS_REPOSITORY_DEFINE_TPLS(
  # Teuchos TPLs
  BinUtils   "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    ST
  ARPREC     "${${PROJECT_NAME}_SOURCE_DIR}/cmake/tpls/"      ST
  QD         "${${PROJECT_NAME}_SOURCE_DIR}/cmake/tpls/"      ST
  MPI        "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"  PT
  BLAS       "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PT
  LAPACK     "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PT
  Boost      "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    ST
  QT         "${${PROJECT_NAME}_SOURCE_DIR}/cmake/tpls/"      ST
  Eigen      "${${PROJECT_NAME}_SOURCE_DIR}/cmake/tpls/"      ST
  # Kokkos TPLs
  Pthread    "${${PROJECT_NAME}_SOURCE_DIR}/kokkos/cmake/tpls/"   ST
  CUDA       "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"      ST
  HWLOC      "${${PROJECT_NAME}_SOURCE_DIR}/kokkos/cmake/tpls/"   ST
  QTHREAD    "${${PROJECT_NAME}_SOURCE_DIR}/kokkos/cmake/tpls/"   ST
  CUSPARSE   "${${PROJECT_NAME}_SOURCE_DIR}/kokkos/cmake/tpls/"   ST
  )
