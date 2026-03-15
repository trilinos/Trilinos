if(NOT CBLAS_ROOT)
  set(CBLAS_ROOT $ENV{CBLAS_ROOT})
endif()
if(NOT CBLAS_ROOT)
  set(CBLAS_ROOT $ENV{OPENBLAS_ROOT})
endif()

if(CBLAS_LIBRARIES)
  #we were given the exact list of libraries to find
  kokkoskernels_find_imported(CBLAS INTERFACE
    LIBRARIES     ${CBLAS_LIBRARIES}
    LIBRARY_PATHS ${CBLAS_LIBRARY_DIRS}
    HEADERS       cblas.h
    HEADER_PATHS  ${CBLAS_INCLUDE_DIRS})
else()
  #we need to find one of the valid versions from the list below
  kokkoskernels_find_imported(CBLAS
    LIBRARY       cblas openblas blas blis
    LIBRARY_PATHS ${CBLAS_LIBRARY_DIRS}
    HEADERS       cblas.h
    HEADER_PATHS  ${CBLAS_INCLUDE_DIRS})
endif()
