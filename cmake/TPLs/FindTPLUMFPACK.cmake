INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( UMFPACK
  REQUIRED_HEADERS umfpack.h amd.h UFconfig.h
  REQUIRED_LIBS_NAMES umfpack amd
  )
