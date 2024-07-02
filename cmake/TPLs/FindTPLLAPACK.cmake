
if (MSVC AND NOT
    (LAPACK_LIBRARY_DIRS  OR
     (NOT "${LAPACK_LIBRARY_NAMES}" STREQUAL "lapack lapack_win32" AND
      NOT "${LAPACK_LIBRARY_NAMES}" STREQUAL "") OR
     LAPACK_INCLUDE_DIRS  OR
     LAPACK_INCLUDE_NAMES OR
     (NOT "${TPL_LAPACK_LIBRARIES}" STREQUAL "lapack" AND
      NOT "${TPL_LAPACK_LIBRARIES}" STREQUAL "") OR
     TPL_LAPACK_INCLUDE_DIRS)
   )
  if(CLAPACK_FOUND)
    advanced_set(TPL_LAPACK_LIBRARIES lapack
        CACHE FILEPATH "Set from MSVC CLAPACK specialization")
  endif()
endif()

tribits_tpl_find_include_dirs_and_libraries( LAPACK
  REQUIRED_LIBS_NAMES "lapack lapack_win32")
