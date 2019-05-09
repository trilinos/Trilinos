#
# cmake -P script that calls tribits_ctest_update_commands.cmake and sends
# STDOUT to a file.
#
# This script is required if you want to capture the output from these
# commands to a file since this is called from ctest_update() which discards
# the output (and does not send it to CDash).
#

message("\ncmake -P tribits_ctest_update_commands_wrapper.cmake:")
message("-- OUTPUT_FILE=${OUTPUT_FILE}\n")

execute_process(
  COMMAND "${CMAKE_COMMAND}"
    -DGIT_EXE=${GIT_EXE}
    -DREMOTE_NAME=${REMOTE_NAME}
    -DBRANCH=${BRANCH}
    -DUNIT_TEST_MODE=${UNIT_TEST_MODE}
    -P ${CMAKE_CURRENT_LIST_DIR}/tribits_ctest_update_commands.cmake
  OUTPUT_FILE "${OUTPUT_FILE}"
  ERROR_FILE "${OUTPUT_FILE}"
  RESULT_VARIABLE RTN_CODE
  )
message("\ntribits_ctest_update_commands_wrapper.cmake return: ${RTN_CODE}\n")

if (NOT "${RTN_CODE}" STREQUAL "0")
  message(FATAL_ERROR "Git Update FAILED!")
endif()
