################################################################################
#
# Common CTest -S script for all test builds on shiller/hansen
#
################################################################################

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../TrilinosCTestDriverCore.atdm.cmake")

# Add this script and the shiller env script to the notes
SET( CTEST_NOTES_FILES
  ${CTEST_NOTES_FILES}
  "${TRIBITS_PROJECT_ROOT}/cmake/std/atdm/shiller/environment.sh"
  "${CMAKE_CURRENT_LIST_FILE}"
  )

# Shiller and Hansen nodes have 16 cores which can support 2 threads each.
# For building, we use all 32 threads.  For running tests, we assume to use
# all of the threads as well with 2 threads per core with an OpenMP build.
SET(CTEST_PARALLEL_LEVEL "16")  # Should be set to 32 if using only 1 thread?

IF ($ENV{ATDM_CONFIG_USE_MAKEFILES})
  SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  SET(CTEST_BUILD_FLAGS "-j32 -k")
ELSE()
  SET(CTEST_CMAKE_GENERATOR Ninja)
  SET(CTEST_BUILD_FLAGS "-k 999999")
ENDIF()

SET(EXTRA_CONFIGURE_OPTIONS)

# Must set the site name so that it does not change depending on what node
# runs the build.
SET(CTEST_SITE "hansen/shiller")

# Run the genetic ATDM driver
TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
