if(NOT LAPACKE_ROOT)
  set(LAPACKE_ROOT $ENV{LAPACKE_ROOT})
endif()
if(NOT LAPACKE_ROOT)
  set(LAPACKE_ROOT $ENV{OPENBLAS_ROOT})
endif()
if(LAPACKE_LIBRARIES)
  #we were given the exact list of libraries to find
  kokkoskernels_find_imported(LAPACKE INTERFACE
    LIBRARIES     ${LAPACKE_LIBRARIES}
    LIBRARY_PATHS ${LAPACKE_LIBRARY_DIRS}
    HEADERS       lapacke.h
    HEADER_PATHS  ${LAPACKE_INCLUDE_DIRS})
else()
  #we need to find one of the valid versions from the list below
  kokkoskernels_find_imported(LAPACKE
    LIBRARY       lapacke openblas
    LIBRARY_PATHS ${LAPACKE_LIBRARY_DIRS}
    HEADERS       lapacke.h
    HEADER_PATHS  ${LAPACKE_INCLUDE_DIRS})
endif()
