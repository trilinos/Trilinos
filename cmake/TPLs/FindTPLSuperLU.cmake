set(REQUIRED_HEADERS supermatrix.h slu_ddefs.h )
set(REQUIRED_LIBS_NAMES "superlu superlu_3.0 superlu_4.0 superlu_4.1 superlu_4.2 superlu_4.3 superlu_5.0 superlu_5.1.1 superlu_5.2.1")
set(IMPORTED_TARGETS_FOR_ALL_LIBS superlu::superlu )

tribits_tpl_allow_pre_find_package(SuperLU SUPERLU_ALLOW_PREFIND)

if (SUPERLU_ALLOW_PREFIND)
  message("-- Using find_package(SuperLU ...) ...")
  find_package(superlu)
  if (superlu_FOUND)
    message("-- Found superlu_DIR='${superlu_DIR}'")
    message("-- Generating SuperLU::all_libs and SuperLUConfig.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(SuperLU
      INNER_FIND_PACKAGE_NAME  superlu
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS} )
  endif()
endif()

if (NOT TARGET SuperLU::all_libs)

  TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( SuperLU
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )

endif()

include(CheckCSourceCompiles)
include(MultilineSet)

# API change in SuperLU 5.0 and 5.1 requires a 'GlobalLU_t' parameter for
# *gssvx, *gsisx, *gstrf, and *gsitrf routines.  Check whether these
# parameters are needed.

FUNCTION(CHECK_SUPERLU_GLOBALLU_T_ARG  VARNAME)
  SET(SOURCE
  "
#include <slu_ddefs.h>

int main()
{
  GlobalLU_t lu;
  superlu_options_t opt;
  SuperMatrix M;
  int *i;
  double *d;
  void *v;
  char *c;
  SuperLUStat_t stat;
  mem_usage_t mem;

  dgsisx(&opt,&M,i,i,i,c,d,d,&M,&M,v,*i,&M,&M,d,d,&lu,&mem,&stat,i);
  return 0;
}
"
  )

  SET(CMAKE_REQUIRED_LIBRARIES SuperLU::all_libs)
  CHECK_C_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

IF (TPL_ENABLE_SuperLU)
  CHECK_SUPERLU_GLOBALLU_T_ARG(${PROJECT_NAME}_ENABLE_SuperLU5_API)
ENDIF(TPL_ENABLE_SuperLU)
