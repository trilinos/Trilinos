#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Belos::belos" for configuration "RELEASE"
set_property(TARGET Belos::belos APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Belos::belos PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libbelos.a"
  )

list(APPEND _cmake_import_check_targets Belos::belos )
list(APPEND _cmake_import_check_files_for_Belos::belos "${_IMPORT_PREFIX}/lib/libbelos.a" )

# Import target "Belos::belosepetra" for configuration "RELEASE"
set_property(TARGET Belos::belosepetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Belos::belosepetra PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libbelosepetra.a"
  )

list(APPEND _cmake_import_check_targets Belos::belosepetra )
list(APPEND _cmake_import_check_files_for_Belos::belosepetra "${_IMPORT_PREFIX}/lib/libbelosepetra.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
