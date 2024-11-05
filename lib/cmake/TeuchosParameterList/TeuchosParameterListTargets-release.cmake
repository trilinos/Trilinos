#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TeuchosParameterList::teuchosparameterlist" for configuration "RELEASE"
set_property(TARGET TeuchosParameterList::teuchosparameterlist APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TeuchosParameterList::teuchosparameterlist PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libteuchosparameterlist.a"
  )

list(APPEND _cmake_import_check_targets TeuchosParameterList::teuchosparameterlist )
list(APPEND _cmake_import_check_files_for_TeuchosParameterList::teuchosparameterlist "${_IMPORT_PREFIX}/lib/libteuchosparameterlist.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
