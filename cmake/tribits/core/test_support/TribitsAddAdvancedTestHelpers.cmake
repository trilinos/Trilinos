# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include("${CMAKE_CURRENT_LIST_DIR}/TribitsAddTestHelpers.cmake")


# Set default ax number of TEST_<idx> blocks in tribits_add_advanced_test()
if ("${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS}" STREQUAL "")
  set(TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS 20)
endif()
# NOTE: Given how includes are done in CMake, above is the only safe way to set
# the default.


# Compute TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX given
# TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS
macro(tribits_add_advanced_test_max_num_test_cmnd_idx_compute)
  math(EXPR TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX
    "${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS}-1")
endmacro()


# Check that the numbers of TEST_<idx> blocks is not violated.
#
# This gets called on the last # TEST_<idx> block of arguments in the parser.
#
function(tribits_add_advanced_test_check_exceed_max_num_test_blocks)
  if (NOT "${PARSE_TEST_${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX}}"
    STREQUAL ""
    )
    list(FIND PARSE_TEST_${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX}
      "TEST_${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS}"
      TEST_BLOCK_KEYWORD_FIND_IDX)
    if (NOT TEST_BLOCK_KEYWORD_FIND_IDX EQUAL -1)
      message_wrapper(FATAL_ERROR
        "${TEST_NAME}: ERROR: Test block TEST_${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS} exceeds the max allowed test block TEST_${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX} as allowed by TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS=${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS}.  To fix this, call set(TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS <larger-num>) before calling tribits_add_advanced_test()!")
      # NOTE: The above error message is put on one line to simplify unit
      # testing checking the string.
    endif()
  endif()
endfunction()


# Join arguments together to add to a set(...) statement to passed to
# execute_process(...)
#
function(tribits_join_exec_process_set_args  OUTPUT_STRING_VAR)
  set(OUTPUT_STRING "")
  foreach(STRING_VAL_RAW ${ARGN})
    # Remove quotes around arguments because CTest does not need them
    string(REGEX REPLACE "\"" "" STRING_VAL "${STRING_VAL_RAW}")
    if (OUTPUT_STRING STREQUAL "")
      set(OUTPUT_STRING "\"${STRING_VAL}\"")
    else()
      set(OUTPUT_STRING "${OUTPUT_STRING} \"${STRING_VAL}\"")
    endif()
  endforeach()
  set(${OUTPUT_STRING_VAR} "${OUTPUT_STRING}" PARENT_SCOPE)
endfunction()


# Unit test helper function for tribits_add_advanced_test(...) that resets
# state before calling tribits_add_advanced_test(...) in unit test mode.
#
# NOTE: The variables:
#
#   TRIBITS_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX}
#
# may get set, even if the test does not get added.  That is because the
# decision to add a test may not be made until much later in the script.
# Therefore, to determine if an advanced test is added or not, one should
# check the variable:
#
#   TRIBITS_ADD_ADVANCED_TEST_NUM_CMNDS
#
# If that variable is empty, then the advanced test did *not* get added.
#
function(tribits_add_advanced_test_unittest_reset)

  global_set(TRIBITS_ADD_ADVANCED_TEST_UNITTEST TRUE)

  global_set(TRIBITS_ADD_ADVANCED_TEST_NUM_CMNDS "")

  tribits_add_advanced_test_max_num_test_cmnd_idx_compute()
  foreach( TEST_CMND_IDX RANGE
      ${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX}
    )
    global_set(TRIBITS_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX} "")
  endforeach()

endfunction()
