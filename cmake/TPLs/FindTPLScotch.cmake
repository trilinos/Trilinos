IF (TPL_ENABLE_MPI)
   TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Scotch
    REQUIRED_HEADERS ptscotch.h
    REQUIRED_LIBS_NAMES ptscotch ptscotcherr scotch scotcherr
   )
ELSE()
   TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Scotch
    REQUIRED_HEADERS scotch.h
    REQUIRED_LIBS_NAMES scotch scotcherr
   )
ENDIF()


# Zoltan2 has a dependency on Scotch 6.0.3

include(CheckCSourceCompiles)
FUNCTION(CHECK_SCOTCH_VERSION_6_0_3  VARNAME)
  SET(SOURCE
  "
  #include <stdio.h>
  #include <stdint.h>
#ifdef TPL_ENABLE_MPI
  #include <mpi.h>
  #include <ptscotch.h>
#else
  #include <scotch.h>
#endif
  int main()
  {
    #if SCOTCH_VERSION > 6 
      return 0;
    #elif SCOTCH_VERSION == 6 && SCOTCH_RELEASE > 0
      return 0;
    #elif SCOTCH_VERSION == 6 && SCOTCH_RELEASE == 0 && SCOTCH_PATCHLEVEL >= 3
      return 0;
    #else
      scotch_version_failure
    #endif
  }
  "
  )
  SET(CMAKE_REQUIRED_LIBRARIES Scotch::all_libs)
  CHECK_C_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

IF(TPL_ENABLE_Scotch)
  CHECK_SCOTCH_VERSION_6_0_3(HAVE_SCOTCH_VERSION_6_0_3)
ENDIF()
