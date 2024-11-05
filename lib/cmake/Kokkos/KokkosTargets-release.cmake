#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Kokkos::kokkoscore" for configuration "RELEASE"
set_property(TARGET Kokkos::kokkoscore APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Kokkos::kokkoscore PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkokkoscore.a"
  )

list(APPEND _cmake_import_check_targets Kokkos::kokkoscore )
list(APPEND _cmake_import_check_files_for_Kokkos::kokkoscore "${_IMPORT_PREFIX}/lib/libkokkoscore.a" )

# Import target "Kokkos::kokkoscontainers" for configuration "RELEASE"
set_property(TARGET Kokkos::kokkoscontainers APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Kokkos::kokkoscontainers PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkokkoscontainers.a"
  )

list(APPEND _cmake_import_check_targets Kokkos::kokkoscontainers )
list(APPEND _cmake_import_check_files_for_Kokkos::kokkoscontainers "${_IMPORT_PREFIX}/lib/libkokkoscontainers.a" )

# Import target "Kokkos::kokkosalgorithms" for configuration "RELEASE"
set_property(TARGET Kokkos::kokkosalgorithms APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Kokkos::kokkosalgorithms PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkokkosalgorithms.a"
  )

list(APPEND _cmake_import_check_targets Kokkos::kokkosalgorithms )
list(APPEND _cmake_import_check_files_for_Kokkos::kokkosalgorithms "${_IMPORT_PREFIX}/lib/libkokkosalgorithms.a" )

# Import target "Kokkos::kokkossimd" for configuration "RELEASE"
set_property(TARGET Kokkos::kokkossimd APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Kokkos::kokkossimd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkokkossimd.a"
  )

list(APPEND _cmake_import_check_targets Kokkos::kokkossimd )
list(APPEND _cmake_import_check_files_for_Kokkos::kokkossimd "${_IMPORT_PREFIX}/lib/libkokkossimd.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
