tribits_package(TrilinosBuildStats)

include("${CMAKE_CURRENT_LIST_DIR}/BuildStatsGatherTarget.cmake")
add_target_gather_build_stats()
# NOTE: We define this build target here after all of the other packages are
# done getting defined so that it will have a dependency on all defined
# library and executable targets but before the end of the tribits_project()
# command so that the TriBITS Ninja Makefile generator will pick this up in
# the documentation!

#
# Create summary report of build stats if build stats are enabled
#

tribits_add_advanced_test( Results
  EXCLUDE_IF_NOT_TRUE ${PROJECT_NAME}_ENABLE_BUILD_STATS
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  TEST_0
    MESSAGE "Gather up the build stats in case the build failed"
    CMND "${CMAKE_CURRENT_LIST_DIR}/gather_build_stats.py"
    ARGS -v
    WORKING_DIRECTORY "${${PROJECT_NAME}_BINARY_DIR}"
    SKIP_CLEAN_WORKING_DIRECTORY  # Critical or you delete the entire working dir!
    ALWAYS_FAIL_ON_NONZERO_RETURN
  TEST_1
    MESSAGE "Sumarize the build stats from the already created build_stats.csv file (CTEST_FULL_OUTPUT)"
    CMND "${${PROJECT_NAME}_SOURCE_DIR}/commonTools/build_stats/summarize_build_stats.py"
    ARGS --bin-by-subdirs-under-dirs=commonTools,packages
      "${${PROJECT_NAME}_BINARY_DIR}/build_stats.csv"
    ALWAYS_FAIL_ON_NONZERO_RETURN
  ADDED_TEST_NAME_OUT Results_TEST_NAME
  )
  # Note, above CTEST_FULL_OUTPUT is in the MESSAGE to get ctest to keep the
  # full test output and send and display on CDash.

if (Results_TEST_NAME)
  set_tests_properties( ${Results_TEST_NAME} PROPERTIES
    ATTACHED_FILES "${${PROJECT_NAME}_BINARY_DIR}/build_stats.csv")
endif()

# NOTE: Above, it is harmless to gather up the build_stats.csv file again in
# this test if the build passed and it has already been gathered up.  But if
# the build failed, then the file build_stats.csv never got updated so it is
# critcial to update it in this test.  Otherwise, you will be displaying an
# old build_stats.csv file from a previous build, which is not good.


#
# Add unit tests for build_stats related support code (even if build stats is
# not enabled for this build).
#

tribits_add_test_directories(unit_tests)


tribits_package_postprocess()
