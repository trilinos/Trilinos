include(CheckCSourceCompiles)

IF (NOT SUPERLU_ROOT)
  SET(SUPERLU_ROOT $ENV{SUPERLU_ROOT})
ENDIF()
IF (SUPERLU_LIBRARIES)
  #we were given the exact list of libraries to find
  KOKKOSKERNELS_FIND_IMPORTED(SUPERLU INTERFACE
    LIBRARIES ${SUPERLU_LIBRARIES}
    LIBRARY_PATHS ${SUPERLU_LIBRARY_DIRS}
    HEADERS slu_ddefs.h
          HEADER_PATHS ${SUPERLU_INCLUDE_DIRS})
ELSE ()
  #we need to find one of the valid versions from the list below
  KOKKOSKERNELS_FIND_IMPORTED(SUPERLU
          LIBRARY superlu
          LIBRARY_PATHS ${SUPERLU_LIBRARY_DIRS}
          HEADERS slu_ddefs.h
          HEADER_PATHS ${SUPERLU_INCLUDE_DIRS})
ENDIF ()

# From Trilinos/cmake/TPLs/FindTPLSuperLU.cmake
FUNCTION(CHECK_SUPERLU_GLOBALLU_T_ARG VARNAME)
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

  SET(CMAKE_REQUIRED_INCLUDES ${TPL_SuperLU_INCLUDE_DIRS})
  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_SuperLU_LIBRARIES} ${TPL_METIS_LIBRARIES} ${TPL_BLAS_LIBRARIES})
  SET(CMAKE_REQUIRED_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
  CHECK_C_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
ENDFUNCTION()

IF (KokkosKernels_ENABLE_TPL_SUPERLU)
  CHECK_SUPERLU_GLOBALLU_T_ARG(${PROJECT_NAME}_ENABLE_SuperLU5_API)
ENDIF (KokkosKernels_ENABLE_TPL_SUPERLU)
