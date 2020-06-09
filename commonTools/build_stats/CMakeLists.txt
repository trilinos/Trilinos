tribits_package(TrilinosBuildStats)

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
    MESSAGE "Sumarize the build stats from the already created build_stats.csv file"
    CMND "${${PROJECT_NAME}_SOURCE_DIR}/commonTools/build_stats/summarize_build_stats.py"
    ARGS --build-stats-csv-file="${${PROJECT_NAME}_BINARY_DIR}/build_stats.csv"
    ALWAYS_FAIL_ON_NONZERO_RETURN
  ADDED_TEST_NAME_OUT Results_TEST_NAME
  )

if (Results_TEST_NAME)
  set_tests_properties( ${Results_TEST_NAME} PROPERTIES
    ATTACHED_FILES "${${PROJECT_NAME}_BINARY_DIR}/build_stats.csv")
endif()


#
# Add unit tests for build_stats related support code (even if build stats is
# not enabled for this build).
#

tribits_add_test_directories(unit_tests)


tribits_package_postprocess()
