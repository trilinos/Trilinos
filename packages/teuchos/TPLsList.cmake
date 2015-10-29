TRIBITS_REPOSITORY_DEFINE_TPLS(
  BinUtils   "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    SS
  ARPREC     "${${PROJECT_NAME}_SOURCE_DIR}/cmake/tpls/"      SS
  QD         "${${PROJECT_NAME}_SOURCE_DIR}/cmake/tpls/"      SS
  MPI        "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"  PS
  BLAS       "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PS
  LAPACK     "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PS
  Boost      "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    SS
  QT         "${${PROJECT_NAME}_SOURCE_DIR}/cmake/tpls/"      ST
  Eigen      "${${PROJECT_NAME}_SOURCE_DIR}/cmake/tpls/"      ST
  )
