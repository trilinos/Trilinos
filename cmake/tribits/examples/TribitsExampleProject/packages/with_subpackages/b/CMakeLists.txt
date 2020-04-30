TRIBITS_SUBPACKAGE(B)

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(tests)

INSTALL( DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/stuff"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/share/${PACKAGE_NAME}"
  USE_SOURCE_PERMISSIONS PATTERN "*~" EXCLUDE )
# Above, we must use 'USE_SOURCE_PERMISSIONS' to preserve the executable
# permission on the scripts in that directory.  (TriBITS will add commands to
# fix the permissions after the install.)

TRIBITS_SUBPACKAGE_POSTPROCESS()
