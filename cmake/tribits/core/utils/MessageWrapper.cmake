# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
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

