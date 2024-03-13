message("+----------------------------------+")
message("| ctest-stage-coverage.cmake START |")
message("+----------------------------------+")
include(${CMAKE_CURRENT_LIST_DIR}/ctest-functions.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)
include(CTestCoverageCollectGCOV)


# -----------------------------------------------------------
# -- Coverage
# -----------------------------------------------------------
set(STAGE_COVERAGE_ERROR OFF)

# Rather than attempt to parse this out via GenConfig and
# pass it through the PR driver scripts, use a substring
# match here.
if(CTEST_BUILD_NAME MATCHES ".*coverage.*")
    banner("START coverage step")

    find_program(GCOV_EXECUTABLE
            NAMES gcov
            DOC "Path to gcov executable")


    # Execute the coverage parsing step
    ctest_coverage_collect_gcov(TARBALL "${CTEST_BUILD_NAME}-gcov.tar.gz"
                                GCOV_COMMAND ${GCOV_EXECUTABLE}
                                GCOV_OPTIONS --branch-probabilities --hash-filenames --demangled-names --human-readable --display-progress
                                GLOB ON)


    #find_program(GCOVR_EXECUTABLE
    #        NAMES gcovr
    #        DOC "Path to gcovr executable")

    # TODO: Try using gcovr
    #execute_process(COMMAND ${GCOVR_EXECUTABLE} --json --output ${CTEST_BUILD_NAME}-gcovr-report.json --exclude-throw-branches --exlude-unreachable-branches -j 29 --gcov-executable ${GCOV_EXECUTABLE} --root ${CTEST_SOURCE_DIRECTORY} ${CTEST_BINARY_DIRECTORY}
    #                COMMAND rm -f ${CTEST_BUILD_NAME}-gcovr-report.tar.gz
    #                COMMAND tar -czf ${CTEST_BUILD_NAME}-gcovr-report.tar.gz ${CTEST_BUILD_NAME}-gcovr-report.json)

    # Print out final stage banner
    banner("END coverage step")

    banner("ctest_upload() START")

    ctest_submit(CDASH_UPLOAD ${CTEST_BUILD_NAME}-gcov.tar.gz
                CDASH_UPLOAD_TYPE GcovTar
                RETRY_COUNT ${ctest_submit_retry_count}
                RETRY_DELAY ${ctest_submit_retry_delay}
                RETURN_VALUE coverage_upload_error)

    if(coverage_upload_error EQUAL 0)
        message(">>> Coverage Upload: OK")
    else()
        message(">>> Coverage Upload: FAILED")
        message(">>> - The ERROR code is ${coverage_upload_error}")
        set(STAGE_COVERAGE_ERROR ON)
    endif()

    banner("ctest_upload() FINISH")

    # TODO: Check for coverage percentage drop?
    #if( NOT (coverage_error EQUAL 0) )
    #    message(WARNING "Coverage has decreased from X to Y.")
    #    set(STAGE_COVERAGE_ERROR ON)
    #else()
    #    message("Coverage passed.")
    #endif()

    message("+-----------------------------------+")
    message("| ctest-stage-coverage.cmake FINISH |")
    message("+-----------------------------------+")
else()
    message("+------------------------------------+")
    message("| ctest-stage-coverage.cmake SKIPPED |")
    message("+------------------------------------+")
endif()
