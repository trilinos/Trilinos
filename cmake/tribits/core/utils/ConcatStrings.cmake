# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include("${CMAKE_CURRENT_LIST_DIR}/PrintVar.cmake")


# @FUNCTION: concat_strings()
#
# Concatenate a set of string arguments.
#
# Usage::
#
#   concat_strings(<outputVar> "<str0>" "<str1>" ...)
#
# On output, ``<outputVar>`` is set to ``"<str0><str1>..."``.  This makes it
# easier to format a long string over multiple CMake source code lines.
#
function(concat_strings OUTPUT_STRING_VAR)
  #message("CONCAT_STRINGS OUTPUT_STRING_VAR: ${OUTPUT_STRING_VAR} {${ARGN}}")
  #print_var(${OUTPUT_STRING_VAR})
  set(OUTPUT_STRING "")
  #print_var(OUTPUT_STRING)
  foreach(STRING_VAL ${ARGN})
    #print_var(STRING_VAL)
    set(OUTPUT_STRING "${OUTPUT_STRING}${STRING_VAL}")
    #print_var(OUTPUT_STRING)
  endforeach()
  set(${OUTPUT_STRING_VAR} "${OUTPUT_STRING}" PARENT_SCOPE)
endfunction()
