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

include("${CMAKE_CURRENT_LIST_DIR}/ConcatStrings.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/PrintVar.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/TribitsDeprecatedHelpers.cmake")


# @FUNCTION: append_string_var()
#
# Append strings to an existing string variable (reduces boiler-place code and
# reduces mistakes).
#
# Usage::
#
#   append_string_var(<stringVar> "<string1>" "<string2>" ...)
#
# Note that the usage of the characters ``'['``, ``']'``, ``'{'``, ``'}'`` are
# taken by CMake to bypass the meaning of ';' to separate string characters.
# If one wants to ignore the meaning of these special characters and are okay
# with just adding one string at a time, then use `append_string_var_ext()`_.
#
# **DEPRECATED**: Instead, use::
#
#   string(APPEND <stringVar> "<string1>" "<string2>" ...)
#
function(append_string_var STRING_VAR_OUT)
  tribits_deprecated_command(append_string_var
    MESSAGE "Use string(APPEND) instead.")
  #message("APPEND_STRING_VAR: ${STRING_VAR_OUT} {${ARGN}}")
  concat_strings( STRING_VAR "${${STRING_VAR_OUT}}" ${ARGN} )
  #print_var( STRING_VAR )
  set(${STRING_VAR_OUT} "${STRING_VAR}" PARENT_SCOPE)
  #print_var(STRING_VAR_OUT)
endfunction()


# @FUNCTION: append_string_var_ext()
#
# Append a single string to an existing string variable, ignoring ';' (reduces
# boiler-place code and reduces mistakes).
#
# Usage::
#
#   append_string_var_ext(<stringVar> "<string>")
#
# Simply sets ``<stringVar> = "${<stringVar>}<string>"`` and leaves in ``';'``
# without creating new array elements.
#
function(append_string_var_ext  STRING_VAR_OUT  STRING_TO_APPEND)
  set(STRING_VAR "${${STRING_VAR_OUT}}${STRING_TO_APPEND}")
  set(${STRING_VAR_OUT} "${STRING_VAR}" PARENT_SCOPE)
endfunction()
