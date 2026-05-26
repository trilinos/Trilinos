if(NOT CHOLMOD_ROOT)
  set(CHOLMOD_ROOT $ENV{SUITESPARSE_ROOT})
endif()
if(CHOLMOD_LIBRARIES)
  #we were given the exact list of libraries to find
  kokkoskernels_find_imported(CHOLMOD INTERFACE
    LIBRARIES     ${CHOLMOD_LIBRARIES}
    LIBRARY_PATHS ${CHOLMOD_LIBRARY_DIRS}
    HEADERS       cholmod.h
    HEADER_PATHS  ${CHOLMOD_INCLUDE_DIRS})
else()
  #we need to find one of the valid versions from the list below
  kokkoskernels_find_imported(CHOLMOD
    LIBRARY       cholmod
    LIBRARY_PATHS ${CHOLMOD_LIBRARY_DIRS}
    HEADERS       cholmod.h
    HEADER_PATHS  ${CHOLMOD_INCLUDE_DIRS})
endif()
