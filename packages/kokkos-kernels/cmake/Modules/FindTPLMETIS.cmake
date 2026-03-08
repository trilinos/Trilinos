if(NOT METIS_ROOT)
  set(METIS_ROOT $ENV{METIS_ROOT})
endif()
#we need to find one of the valid versions from the list below
kokkoskernels_find_imported(METIS
  LIBRARY       metis
  LIBRARY_PATHS ${METIS_LIBRARY_DIRS}
  HEADERS       metis.h
  HEADER_PATHS  ${METIS_INCLUDE_DIRS})
