tribits_subpackage(B)

add_subdirectory(src)

# Set variables that will be used in this package and downstream packages

if (${PROJECT_NAME}_ENABLE_InsertedPkg)
  global_set(EXPECTED_INSERTEDPKG_AND_DEPS "InsertedPkg ${EXPECTED_SIMPLECXX_AND_DEPS}")
  global_set(EXPECTED_INSERTEDPKG_AND_DEPS_STR "${EXPECTED_INSERTEDPKG_AND_DEPS} ")
else()
  global_set(EXPECTED_INSERTEDPKG_DEPS "")
  global_set(EXPECTED_INSERTEDPKG_DEPS_STR "")
endif()

if (${PROJECT_NAME}_ENABLE_WithSubpackagesA)
  global_set(EXPECTED_A_AND_DEPS "A ${EXPECTED_SIMPLECXX_AND_DEPS}")
  global_set(EXPECTED_A_AND_DEPS_STR "${EXPECTED_A_AND_DEPS} ")
else()
  global_set(EXPECTED_A_AND_DEPS "")
  global_set(EXPECTED_A_AND_DEPS_STR "")
endif()

global_set(EXPECTED_B_DEPS
  "${EXPECTED_A_AND_DEPS_STR}${EXPECTED_INSERTEDPKG_AND_DEPS_STR}${EXPECTED_SIMPLECXX_AND_DEPS}")
tribits_pkg_export_cache_var(EXPECTED_B_DEPS)

tribits_add_test_directories(tests)

install( DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/stuff"
  DESTINATION "share/${PACKAGE_NAME}"
  USE_SOURCE_PERMISSIONS PATTERN "*~" EXCLUDE )
# Above, we must use 'USE_SOURCE_PERMISSIONS' to preserve the executable
# permission on the scripts in that directory.  (TriBITS will add commands to
# fix the permissions after the install.)

tribits_subpackage_postprocess()
