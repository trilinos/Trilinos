message("+--------------------------------------+")
message("| ctest-driver.cmake START             |")
message("+--------------------------------------+")
message("Snapshot: 2021-09-22 001")

include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/ctest-stage-configure.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/ctest-stage-build.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/ctest-stage-test.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/ctest-stage-coverage.cmake)

if(STAGE_CONFIGURE_ERROR OR STAGE_BUILD_ERROR OR STAGE_TEST_ERROR OR STAGE_COVERAGE_ERROR)
    message(FATAL_ERROR "STAGE_CONFIGURE_ERROR: ${STAGE_CONFIGURE_ERROR}, STAGE_BUILD_ERROR: ${STAGE_BUILD_ERROR}, STAGE_TEST_ERROR: ${STAGE_TEST_ERROR}, STAGE_COVERAGE_ERROR: ${STAGE_COVERAGE_ERROR}")
endif()

message(">>> CDash URL1 = ${build_url1}")
message(">>> CDash URL2 = ${build_url2}")
message(">>> CDash URL3 = ${build_url3}")

message("+--------------------------------------+")
message("| ctest-driver.cmake FINISH            |")
message("+--------------------------------------+")
