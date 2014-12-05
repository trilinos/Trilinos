TRIBITS_REPOSITORY_DEFINE_TPLS(
  MPI            "${${PROJECT_NAME}_TRIBITS_DIR}/tpls/"    PT
  BLAS           cmake/TPLs/    PT
  LAPACK         cmake/TPLs/    PT
  Boost          cmake/TPLs/    ST
  UMFPACK        cmake/TPLs/    ST
  AMD            cmake/TPLs/    EX
  PETSC          cmake/TPLs/    ST
  )
