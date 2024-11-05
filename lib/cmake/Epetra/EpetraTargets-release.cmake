#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Epetra::epetra" for configuration "RELEASE"
set_property(TARGET Epetra::epetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Epetra::epetra PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX;Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libepetra.a"
  )

list(APPEND _cmake_import_check_targets Epetra::epetra )
list(APPEND _cmake_import_check_files_for_Epetra::epetra "${_IMPORT_PREFIX}/lib/libepetra.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
