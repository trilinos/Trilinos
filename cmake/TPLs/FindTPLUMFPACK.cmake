
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( UMFPACK
  REQUIRED_HEADERS umfpack.h amd.h
  REQUIRED_LIBS_NAMES umfpack amd
  )

# Amesos2 has Umfpack wrappers which depend on being able to
# send complex as Real-Im as a single array.
# Old versions of Umfpack don't support this so for now we are
# just disabling the complex tests for Umfpack version < 5.

include(CheckCSourceCompiles)
FUNCTION(CHECK_UMFPACK_HAS_VERSION_5  VARNAME)
  SET(SOURCE
  "
  #include <stdio.h>
  #include <umfpack.h>
  int main()
  {
    // this got added in Umfpack version 4.5
    #if UMFPACK_MAIN_VERSION >= 5
      return 0;
    #else
      umfpack_version_failure
    #endif
  }

  "
  )
  SET(CMAKE_REQUIRED_INCLUDES ${TPL_UMFPACK_INCLUDE_DIRS})
  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_UMFPACK_LIBRARIES})
  SET(CMAKE_REQUIRED_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
  CHECK_C_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

IF(TPL_ENABLE_UMFPACK)
  CHECK_UMFPACK_HAS_VERSION_5(HAVE_UMFPACK_VERSION_5)
ENDIF()
