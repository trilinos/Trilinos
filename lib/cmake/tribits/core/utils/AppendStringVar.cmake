# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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
