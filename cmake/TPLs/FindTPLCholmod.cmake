
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Cholmod
  REQUIRED_HEADERS cholmod.h cholmod_core.h
  REQUIRED_LIBS_NAMES libcholmod.a libamd.a libcolamd.a libccolamd.a libcamd.a libsuitesparseconfig.a
  )


