#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TeuchosRemainder::teuchosremainder" for configuration "RELEASE"
set_property(TARGET TeuchosRemainder::teuchosremainder APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TeuchosRemainder::teuchosremainder PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libteuchosremainder.a"
  )

list(APPEND _cmake_import_check_targets TeuchosRemainder::teuchosremainder )
list(APPEND _cmake_import_check_files_for_TeuchosRemainder::teuchosremainder "${_IMPORT_PREFIX}/lib/libteuchosremainder.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
