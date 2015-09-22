TRIBITS_REPOSITORY_DEFINE_TPLS(
  MPI        "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"  PS
  BLAS       "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PS
  LAPACK     "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"    PS
  METIS      "dummy_path/"   EX
  ParMETIS   "dummy_path/"   EX
  PETSC      "dummy_path/"   EX
  SuperLU    "dummy_path/"   EX
  MATLAB     "dummy_path/"   EX
  )
