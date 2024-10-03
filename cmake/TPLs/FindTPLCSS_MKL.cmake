
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( CSS_MKL
  REQUIRED_HEADERS mkl_pardiso.h mkl_cluster_sparse_solver.h
  REQUIRED_LIBS_NAMES mkl_intel_lp64
  )
