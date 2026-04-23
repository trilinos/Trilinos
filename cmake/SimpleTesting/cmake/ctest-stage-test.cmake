message("+--------------------------------------+")
message("| ctest-stage-test.cmake START         |")
message("+--------------------------------------+")
include(${CMAKE_CURRENT_LIST_DIR}/ctest-functions.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)

macro(trilinos_ctest_compute_asan_memcheck_route)
  set(TRILINOS_CTEST_USE_ASAN_MEMCHECK FALSE)
  set(TRILINOS_CTEST_USE_ASAN_MEMCHECK_REASON "")
  if(DEFINED CTEST_BUILD_NAME AND CTEST_BUILD_NAME MATCHES ".*_asan_.*")
    set(TRILINOS_CTEST_USE_ASAN_MEMCHECK TRUE)
    set(TRILINOS_CTEST_USE_ASAN_MEMCHECK_REASON "CTEST_BUILD_NAME matches .*_asan_.*")
  endif()
  if(DEFINED CTEST_BINARY_DIRECTORY)
    set(_trilinos_genconfig_bn_file "${CTEST_BINARY_DIRECTORY}/genconfig_build_name.txt")
    if(EXISTS "${_trilinos_genconfig_bn_file}")
      file(READ "${_trilinos_genconfig_bn_file}" _trilinos_genconfig_bn_content)
      string(STRIP "${_trilinos_genconfig_bn_content}" _trilinos_genconfig_bn_content)
      if(_trilinos_genconfig_bn_content MATCHES ".*_asan_.*")
        set(TRILINOS_CTEST_USE_ASAN_MEMCHECK TRUE)
        if(TRILINOS_CTEST_USE_ASAN_MEMCHECK_REASON STREQUAL "")
          set(TRILINOS_CTEST_USE_ASAN_MEMCHECK_REASON
              "genconfig_build_name.txt matches .*_asan_.*")
        else()
          set(TRILINOS_CTEST_USE_ASAN_MEMCHECK_REASON
              "${TRILINOS_CTEST_USE_ASAN_MEMCHECK_REASON}; genconfig_build_name.txt matches .*_asan_.*")
        endif()
      endif()
    endif()
  endif()
endmacro()


# -----------------------------------------------------------
# -- Test
# -----------------------------------------------------------
banner("START test step")

set(STAGE_TEST_ERROR OFF)

trilinos_ctest_compute_asan_memcheck_route()

if(NOT SKIP_RUN_TESTS)
    if(TRILINOS_CTEST_USE_ASAN_MEMCHECK)
        message(">>> AddressSanitizer memcheck route: ${TRILINOS_CTEST_USE_ASAN_MEMCHECK_REASON}")
        set(CTEST_MEMORYCHECK_TYPE "AddressSanitizer")
        set(ENV{LSAN_OPTIONS} "suppressions=${CTEST_SOURCE_DIRECTORY}/packages/framework/asan_assets/lsan.supp")
        set(ENV{LD_PRELOAD} ${CTEST_SOURCE_DIRECTORY}/packages/framework/asan_assets/dummy_dlclose.so)
        ctest_memcheck(PARALLEL_LEVEL ${TEST_PARALLEL_LEVEL}
                       CAPTURE_CMAKE_ERROR captured_cmake_error
                       RETURN_VALUE test_error)
        unset(ENV{LD_PRELOAD})
        submit_by_parts( "MemCheck" )
    else()
        message(">>> Plain ctest_test route (not AddressSanitizer memcheck)")
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

