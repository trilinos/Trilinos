# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()

include(MessageWrapper)

# Called to report incorrect usage
#
# Usage:
#
#   tribits_report_invalid_tribits_usage("<err_msg0>" "<err_msg1>" ...)
#
# Depending on the value of ${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE, this
# function will:
#
#   * FATAL_ERROR: Calls message(FATAL_ERROR "<error_message>")
#   * SEND_ERROR: Calls message(SEND_ERROR "<error_message>")
#   * WARNING: Calls message(WARNING "<error_message>")
#   * NOTICE: Calls message(NOTICE "<error_message>")
#   * IGNORE: Does not call message() at all and is silent
#
# If '${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE}' is empty on call, then
# `FATAL_ERROR` will be used.
#
function(tribits_report_invalid_tribits_usage)
  if ("${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE}" STREQUAL "")
    set(${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE  FATAL_ERROR)
  endif()
  set(printErrMsgMode ${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE})
  set(ignoreValues  "IGNORE" "OFF")
  if(${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE  IN_LIST  ignoreValues)
    set(printErrMsgMode "")
  endif()
  if (printErrMsgMode)
    message_wrapper(${printErrMsgMode} ${ARGN})
  endif()
endfunction()
