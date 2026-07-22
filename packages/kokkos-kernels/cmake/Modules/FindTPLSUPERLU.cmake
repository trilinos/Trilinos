include(CheckCSourceCompiles)

if(NOT SUPERLU_ROOT)
  set(SUPERLU_ROOT $ENV{SUPERLU_ROOT})
endif()
if(SUPERLU_LIBRARIES)
  #we were given the exact list of libraries to find
  kokkoskernels_find_imported(SUPERLU INTERFACE
    LIBRARIES      ${SUPERLU_LIBRARIES}
    LIBRARY_PATHS  ${SUPERLU_LIBRARY_DIRS}
    HEADERS        slu_ddefs.h
    HEADER_PATHS   ${SUPERLU_INCLUDE_DIRS})
else()
  #we need to find one of the valid versions from the list below
  kokkoskernels_find_imported(SUPERLU
    LIBRARY       superlu
    LIBRARY_PATHS ${SUPERLU_LIBRARY_DIRS}
    HEADERS       slu_ddefs.h
    HEADER_PATHS  ${SUPERLU_INCLUDE_DIRS})
endif()

# From Trilinos/cmake/TPLs/FindTPLSuperLU.cmake
function(check_superlu_globallu_t_arg VARNAME)
  set(SOURCE
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
")

  set(CMAKE_REQUIRED_INCLUDES ${TPL_SuperLU_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${TPL_SuperLU_LIBRARIES} ${TPL_METIS_LIBRARIES} ${TPL_BLAS_LIBRARIES})
  set(CMAKE_REQUIRED_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
  check_c_source_compiles("${SOURCE}" ${VARNAME})
endfunction()

if(KokkosKernels_ENABLE_TPL_SUPERLU)
  check_superlu_globallu_t_arg(${PROJECT_NAME}_ENABLE_SuperLU5_API)
endif(KokkosKernels_ENABLE_TPL_SUPERLU)
