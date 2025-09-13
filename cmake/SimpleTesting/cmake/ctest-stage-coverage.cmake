message("+----------------------------------+")
message("| ctest-stage-coverage.cmake START |")
message("+----------------------------------+")
include(${CMAKE_CURRENT_LIST_DIR}/ctest-functions.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)

# -----------------------------------------------------------
# -- Coverage
# -----------------------------------------------------------
set(STAGE_COVERAGE_ERROR OFF)

# Rather than attempt to parse this out via GenConfig and
# pass it through the PR driver scripts, use a substring
# match here.
if(CTEST_BUILD_NAME MATCHES ".*coverage.*")
    banner("START coverage collection step")

    if(NOT DEFINED CTEST_COVERAGE_COMMAND)
        find_program(GCOV_EXECUTABLE
            NAMES gcov
            DOC "Path to gcov executable"
            REQUIRED)
        set(CTEST_COVERAGE_COMMAND "${GCOV_EXECUTABLE}")
    endif()
    message(">>> CTEST_COVERAGE_COMMAND=${CTEST_COVERAGE_COMMAND}")

    ctest_coverage(
        BUILD "${CTEST_BINARY_DIRECTORY}"
        RETURN_VALUE coverage_collect_error)

    if(coverage_collect_error EQUAL 0)
        banner("END coverage collection step - SUCCESS")
    else()
        message(WARNING "Coverage collection step failed with error ${coverage_collect_error}")
        banner("END coverage collection - FAILURE")
        set(STAGE_COVERAGE_ERROR ON)
    endif()

    banner("END coverage step")
    submit_by_parts("Coverage")

    message("+-----------------------------------+")
    message("| ctest-stage-coverage.cmake FINISH |")
    message("+-----------------------------------+")
else()
    message("+------------------------------------+")
    message("| ctest-stage-coverage.cmake SKIPPED |")
    message("+------------------------------------+")
endif()
