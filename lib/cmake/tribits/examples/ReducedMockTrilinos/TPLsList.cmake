tribits_repository_define_tpls(
  MPI            "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"  PT
  BLAS           "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PT
  LAPACK         "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PT
  Boost          cmake/TPLs/    ST
  UMFPACK        cmake/TPLs/    ST
  AMD            cmake/TPLs/    EX
  PETSC          "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    ST
  )
