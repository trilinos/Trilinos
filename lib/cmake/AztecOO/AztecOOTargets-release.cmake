#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "AztecOO::aztecoo" for configuration "RELEASE"
set_property(TARGET AztecOO::aztecoo APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(AztecOO::aztecoo PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX;Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libaztecoo.a"
  )

list(APPEND _cmake_import_check_targets AztecOO::aztecoo )
list(APPEND _cmake_import_check_files_for_AztecOO::aztecoo "${_IMPORT_PREFIX}/lib/libaztecoo.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
