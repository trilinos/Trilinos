message("+--------------------------------------+")
message("| ctest-stage-test.cmake START         |")
message("+--------------------------------------+")
include(${CMAKE_CURRENT_LIST_DIR}/ctest-functions.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)


# -----------------------------------------------------------
# -- Test
# -----------------------------------------------------------
banner("START test step")

set(STAGE_TEST_ERROR OFF)

if(NOT SKIP_RUN_TESTS)
    if(CTEST_BUILD_NAME MATCHES .*_asan_.*)
        set(CTEST_MEMORYCHECK_TYPE "AddressSanitizer")
        set(ENV{LSAN_OPTIONS} "suppressions=${CTEST_SOURCE_DIRECTORY}/packages/framework/asan_assets/lsan.supp")
        set(ENV{LD_PRELOAD} ${CTEST_SOURCE_DIRECTORY}/packages/framework/asan_assets/dummy_dlclose.so)
        ctest_memcheck(PARALLEL_LEVEL ${TEST_PARALLEL_LEVEL}
                       CAPTURE_CMAKE_ERROR captured_cmake_error
                       RETURN_VALUE test_error)
        unset(ENV{LD_PRELOAD})
        submit_by_parts( "MemCheck" )
    else()
        ctest_test(PARALLEL_LEVEL ${TEST_PARALLEL_LEVEL}
                   CAPTURE_CMAKE_ERROR captured_cmake_error
                   RETURN_VALUE test_error)
        submit_by_parts( "Test" )
    endif()
else()
    message(">>> SKIPPED RUNNING TESTS (SKIP_RUN_TESTS=${SKIP_RUN_TESTS})")
    set(test_error 0)
    submit_by_parts("Test")
endif()

# Print out final stage banner
if(${test_error} EQUAL 0)
    banner("END test step - SUCCESS")
else()
    banner("END test step - FAILURE")
endif()

# Upload configure files
submit_upload_config_files()

# Display what return values were captured.
message(">>> captured_cmake_error = ${captured_cmake_error} (unused)")
message(">>> test_error           = ${test_error}")

# We fail the test here if the `test_error` value from ctest_test is
# nonzero.
# We should NOT check `captured_cmake_error` (CAPTURE_CMAKE_ERROR)
# to determine if our test has failed since that value contain a
# nonzero value if there were _no tests to run_. For Trilinos, this
# can happen since we dynamically enable tests based on changes in
# a Pull-Request.
if( NOT (test_error EQUAL 0) )
    message(WARNING "There are errors detected in the test.")
    set(STAGE_TEST_ERROR ON)
else()
    message("Test(s) passed.")
endif()

message("+--------------------------------------+")
message("| ctest-stage-test.cmake FINISH        |")
message("+--------------------------------------+")

