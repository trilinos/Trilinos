# Tertiary stable code since amd is also built as part of UMFPACK.


TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( AMD
  REQUIRED_HEADERS amd.h
  REQUIRED_LIBS_NAMES amd
  )
