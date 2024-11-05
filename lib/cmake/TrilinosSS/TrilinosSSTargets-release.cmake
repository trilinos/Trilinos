#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TrilinosSS::trilinosss" for configuration "RELEASE"
set_property(TARGET TrilinosSS::trilinosss APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TrilinosSS::trilinosss PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtrilinosss.a"
  )

list(APPEND _cmake_import_check_targets TrilinosSS::trilinosss )
list(APPEND _cmake_import_check_files_for_TrilinosSS::trilinosss "${_IMPORT_PREFIX}/lib/libtrilinosss.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
