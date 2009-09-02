# Tertiary stable code since amd is also built as part of UMFPACK.

INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( AMD
  REQUIRED_HEADERS amd.h
  REQUIRED_LIBS_NAMES amd
  )
