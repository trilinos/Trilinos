#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ML::ml" for configuration "RELEASE"
set_property(TARGET ML::ml APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ML::ml PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libml.a"
  )

list(APPEND _cmake_import_check_targets ML::ml )
list(APPEND _cmake_import_check_files_for_ML::ml "${_IMPORT_PREFIX}/lib/libml.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
