#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TeuchosCore::teuchoscore" for configuration "RELEASE"
set_property(TARGET TeuchosCore::teuchoscore APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TeuchosCore::teuchoscore PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libteuchoscore.a"
  )

list(APPEND _cmake_import_check_targets TeuchosCore::teuchoscore )
list(APPEND _cmake_import_check_files_for_TeuchosCore::teuchoscore "${_IMPORT_PREFIX}/lib/libteuchoscore.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
