################################################################################
#
# Common CTest -S script for all ATDM builds
#
################################################################################

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.atdm.cmake")

ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_KNOWN_SYSTEM_NAME)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_USE_MAKEFILES)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_CTEST_PARALLEL_LEVEL)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_BUILD_COUNT)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_SYSTEM_CDASH_SITE)

# Add this script and the shiller env script to the notes
SET( CTEST_NOTES_FILES
  ${CTEST_NOTES_FILES}
  "${TRIBITS_PROJECT_ROOT}/cmake/std/atdm/$ENV{ATDM_CONFIG_KNOWN_SYSTEM_NAME}/environment.sh"
  "${CMAKE_CURRENT_LIST_FILE}"
  )

SET(CTEST_PARALLEL_LEVEL "$ENV{ATDM_CONFIG_CTEST_PARALLEL_LEVEL}")

IF ($ENV{ATDM_CONFIG_USE_MAKEFILES})
  SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  SET(CTEST_BUILD_FLAGS "-j$ENV{ATDM_CONFIG_BUILD_COUNT} -k")
ELSE()
  SET(CTEST_CMAKE_GENERATOR Ninja)
  SET(CTEST_BUILD_FLAGS "-k 999999")
ENDIF()

SET(EXTRA_CONFIGURE_OPTIONS)

# Must set the site name so that it does not change depending on what node
# runs the build.
SET(CTEST_SITE "$ENV{ATDM_CONFIG_SYSTEM_CDASH_SITE}")

# Run the genetic ATDM driver
TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
