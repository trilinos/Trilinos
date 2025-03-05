include(CheckCSourceCompiles)

FUNCTION(ZOLTAN2_CHECK_PARMETIS_HAS_VERSION_4_0_3  VARNAME)
  SET(SOURCE
  "
  #include <stdio.h>
  #include <parmetis.h>
  int main()
  {
    #if PARMETIS_MAJOR_VERSION > 4
      return 0;
    #elif PARMETIS_MAJOR_VERSION == 4 && PARMETIS_MINOR_VERSION > 0
      return 0;
    #elif PARMETIS_MAJOR_VERSION == 4 && PARMETIS_MINOR_VERSION == 0 && PARMETIS_SUBMINOR_VERSION >= 3
      return 0;
    #else
      parmetis_version_failure
    #endif
  }

  "
  )
  SET(CMAKE_REQUIRED_LIBRARIES ParMETIS::all_libs MPI::all_libs)
  CHECK_C_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()
