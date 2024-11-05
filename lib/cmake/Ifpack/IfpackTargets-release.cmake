#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Ifpack::ifpack" for configuration "RELEASE"
set_property(TARGET Ifpack::ifpack APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Ifpack::ifpack PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libifpack.a"
  )

list(APPEND _cmake_import_check_targets Ifpack::ifpack )
list(APPEND _cmake_import_check_files_for_Ifpack::ifpack "${_IMPORT_PREFIX}/lib/libifpack.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
