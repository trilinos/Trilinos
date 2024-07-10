
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( SuperLUMT
  REQUIRED_HEADERS supermatrix.h slu_mt_util.h pssp_defs.h pdsp_defs.h pcsp_defs.h pzsp_defs.h
  REQUIRED_LIBS_NAMES "superlumt superlu_mt superlu_mt_PTHREAD superlu_mt_DEC superlu_mt_altix superlu_mt_CRAY superlu_mt_sp superlu_mt_ORIGIN superlu_mt_SGI superlu_mt_SOLARIS"
  )
