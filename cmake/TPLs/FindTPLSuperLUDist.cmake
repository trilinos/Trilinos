INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( SuperLUDist
  REQUIRED_HEADERS "superlu_defs.h superludefs.h" supermatrix.h
  REQUIRED_LIBS_NAMES "superludist superlu_dist superlu_dist_2.0 superlu_dist_2.5"
  )
