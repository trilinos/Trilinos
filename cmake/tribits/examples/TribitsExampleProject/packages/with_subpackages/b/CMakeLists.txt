tribits_subpackage(B)

add_subdirectory(src)

tribits_add_test_directories(tests)

install( DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/stuff"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/share/${PACKAGE_NAME}"
  USE_SOURCE_PERMISSIONS PATTERN "*~" EXCLUDE )
# Above, we must use 'USE_SOURCE_PERMISSIONS' to preserve the executable
# permission on the scripts in that directory.  (TriBITS will add commands to
# fix the permissions after the install.)

tribits_subpackage_postprocess()
