TRIBITS_REPOSITORY_DEFINE_TPLS(
  MPI            "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"    PT
  BLAS           "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"    PT
  LAPACK         "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"    PT
  Boost          "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"    ST
  UMFPACK        cmake/TPLs/    ST
  AMD            cmake/TPLs/    EX
  PETSC          "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    ST
  )
