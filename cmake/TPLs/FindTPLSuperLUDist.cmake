
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( SuperLUDist
  REQUIRED_HEADERS "superlu_defs.h superludefs.h" supermatrix.h
  REQUIRED_LIBS_NAMES "superludist superlu_dist superlu_dist_2.0 superlu_dist_2.5 superlu_dist_4.0"
  )

include(CheckCSourceCompiles)
include(MultilineSet)

# Versions 3.0 and later of SuperLU_DIST namespace the IterRefine_t
# enum values with "SLU_" (e.g. "SLU_DOUBLE" versus "DOUBLE" in
# previous versions).  Check which style is used so that packages like
# Amesos and Amesos2 can use the correct enum value.
FUNCTION(CHECK_SUPERLUDIST_ENUM_NAMESPACE  VARNAME)
  SET(SOURCE
  "
#include <superlu_enum_consts.h>

int main()
{
  IterRefine_t refine = SLU_DOUBLE;
  return 0;
}
"
  )

  SET(CMAKE_REQUIRED_LIBRARIES SuperLUDist::all_libs)
  SET(CMAKE_REQUIRED_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

# Version 4.0 of SuperLU_DIST changed the calling parameters of the
# LUstructInit function.  Check which is used here.
FUNCTION(CHECK_SUPERLUDIST_LUSTRUCTINIT  VARNAME)
  SET(SOURCE
  "
#include <superlu_ddefs.h>

int main()
{
  LUstruct_t lu;
  /* This will fail to compile if the 3-arg version is declared. */
  LUstructInit(10, &lu);
}
"
  )

  SET(CMAKE_REQUIRED_LIBRARIES SuperLUDist::all_libs)
  SET(CMAKE_REQUIRED_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

IF (TPL_ENABLE_SuperLUDist)
  CHECK_SUPERLUDIST_ENUM_NAMESPACE(HAVE_SUPERLUDIST_ENUM_NAMESPACE)
  CHECK_SUPERLUDIST_LUSTRUCTINIT(HAVE_SUPERLUDIST_LUSTRUCTINIT_2ARG)
ENDIF(TPL_ENABLE_SuperLUDist)
