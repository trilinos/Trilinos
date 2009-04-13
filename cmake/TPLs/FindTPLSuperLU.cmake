INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( SuperLU
  REQUIRED_HEADERS supermatrix.h
  REQUIRED_LIBS_NAMES "superlu superlu_3.0"
  )
