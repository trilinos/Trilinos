# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include("${CMAKE_CURRENT_LIST_DIR}/GlobalSet.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/AssertDefined.cmake")


# @FUNCTION: append_global_set()
#
# Utility macro that appends arguments to a global variable (reduces
# boiler-plate code and mistakes).
#
# Usage::
#
#   append_global_set(<varName> <arg0> <arg1> ...)
#
# NOTE: The variable ``<varName>`` must exist before calling this function.
# To set it empty initially use `global_null_set()`_.
#
function(append_global_set  VARNAME)
  assert_defined(${VARNAME})
  list(APPEND ${VARNAME} ${ARGN})
  global_set(${VARNAME} ${${VARNAME}})
endfunction()
