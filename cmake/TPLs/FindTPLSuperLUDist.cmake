INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( SuperLUDist
  REQUIRED_HEADERS superludefs.h supermatrix.h
  REQUIRED_LIBS_NAMES "superludist superlu_dist superlu_dist_2.0"
  )
