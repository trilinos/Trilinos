INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( Scotch
  REQUIRED_HEADERS scotch.h ptscotch.h
  REQUIRED_LIBS_NAMES ptscotch ptscotcherr ptscotcherrexit
  )
