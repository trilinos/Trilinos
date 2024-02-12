# A utility driver for submitting partial results to CDash
# when an exception of failure occurs.
message("+--------------------------------------+")
message("| ctest-submit.cmake START             |")
message("+--------------------------------------+")
message("Snapshot: 2021-09-22 001")

include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)

message(" -- Submit - ${CTEST_BUILD_NAME} --")

ctest_start(APPEND)

ctest_submit(RETURN_VALUE res)

message("+--------------------------------------+")
message("| ctest-submit.cmake FINISH            |")
message("+--------------------------------------+")
