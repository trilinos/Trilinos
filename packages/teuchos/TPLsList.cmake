TRIBITS_REPOSITORY_DEFINE_TPLS(
  BinUtils   "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    SS
  ARPREC     "dummy_path/"    SS
  QD         "dummy_path/"    SS
  MPI        "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"  PS
  BLAS       "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PS
  LAPACK     "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PS
  Boost      "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    SS
  QT         "dummy_path/"    SS
  Eigen      "dummy_path/"  EX  # Not in TriBITS common_tpls/
  )
