INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( HIPS
  REQUIRED_HEADERS hips.h hips_fortran.h
  REQUIRED_LIBS_NAMES hips hipssequential io spkit
  )
