# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()

include("${CMAKE_CURRENT_LIST_DIR}/GlobalSet.cmake")

# @FUNCTION: message_wrapper()
#
# Function that wraps the standard CMake/CTest ``message()`` function call in
# order to allow unit testing to intercept the output.
#
# Usage::
#
#   message_wrapper(...)
#
# This function takes exactly the same arguments as built-in ``message()``
# function.  However, when the variable ``MESSAGE_WRAPPER_UNIT_TEST_MODE`` is
# set to ``TRUE``, then this function will not call ``message(...)`` but
# instead will prepend set to the global variable ``MESSAGE_WRAPPER_INPUT``
# the input argument that would have gone to ``message()``.  To capture just
# this call's input, first call::
#
#   global_null_set(MESSAGE_WRAPPER_INPUT)
#
# before calling this function (or the functions/macros that call this
# function).
#
# This function allows one to unit test other user-defined CMake macros and
# functions that call this function to catch error conditions without stopping
# the CMake program.  Otherwise, this is used to capture print messages to
# verify that they say the right thing.
#
function(message_wrapper)
  cmake_parse_arguments(PARSE_ARGV 0 FWD "" "" "")
  #message("MESSAGE_WRAPPER: ${FWD_UNPARSED_ARGUMENTS}")
  if (MESSAGE_WRAPPER_UNIT_TEST_MODE)
    global_set(MESSAGE_WRAPPER_INPUT "${MESSAGE_WRAPPER_INPUT}"
      ${FWD_UNPARSED_ARGUMENTS})
  else()
    message(${FWD_UNPARSED_ARGUMENTS})
  endif()
endfunction()

