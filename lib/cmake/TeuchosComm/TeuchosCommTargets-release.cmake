#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TeuchosComm::teuchoscomm" for configuration "RELEASE"
set_property(TARGET TeuchosComm::teuchoscomm APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TeuchosComm::teuchoscomm PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libteuchoscomm.a"
  )

list(APPEND _cmake_import_check_targets TeuchosComm::teuchoscomm )
list(APPEND _cmake_import_check_files_for_TeuchosComm::teuchoscomm "${_IMPORT_PREFIX}/lib/libteuchoscomm.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
