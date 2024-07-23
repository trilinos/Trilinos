# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include("${CMAKE_CURRENT_LIST_DIR}/ConcatStrings.cmake")


# @FUNCTION: append_string_var_with_sep()
#
# Append strings to a given string variable, joining them using a separator
# string.
#
# Usage::
#
#   append_string_var_with_sep(<stringVar> "<sepStr>" "<str0>" "<str1>" ...)
#
# Each of the strings ``<stri>`` are appended to ``<stringVar>`` using the
# separation string ``<sepStr>``.
#
function(append_string_var_with_sep  STRING_VAR  SEP_STR)
  #message("APPEND_STRING_VAR: '${STRING_VAR}' '${SEP_STR}' ${ARGN}")
  #print_var(STRING_VAR)
  #print_var(${STRING_VAR})
  if (${STRING_VAR})
    concat_strings( TMP_STRING "${${STRING_VAR}}${SEP_STR}" ${ARGN} )
  else()
    concat_strings( TMP_STRING ${ARGN} )
  endif()
  #print_var( TMP_STRING )
  set(${STRING_VAR} "${TMP_STRING}" PARENT_SCOPE)
  #set(${STRING_VAR} "${${STRING_VAR}}${LINE}" PARENT_SCOPE)
  #print_var(STRING_VAR)
endfunction()
