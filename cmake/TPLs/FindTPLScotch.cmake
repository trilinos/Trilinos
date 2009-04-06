INCLUDE(TPLDeclareLibraries)

IF (TPL_ENABLE_MPI)
   TPL_DECLARE_LIBRARIES( Scotch
    REQUIRED_HEADERS ptscotch.h
    REQUIRED_LIBS_NAMES ptscotch ptscotcherr
   )
ELSE()
   TPL_DECLARE_LIBRARIES( Scotch
    REQUIRED_HEADERS scotch.h
    REQUIRED_LIBS_NAMES scotch scotcherr
   )
ENDIF()
