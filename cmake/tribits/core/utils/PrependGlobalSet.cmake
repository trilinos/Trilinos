# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(GlobalSet)
include(AssertDefined)


# @MACRO: prepend_global_set()
#
# Utility macro that prepends arguments to a global variable (reduces
# boiler-plate code and mistakes).
#
# Usage::
#
#   prepend_global_set(<varName> <arg0> <arg1> ...)
#
# The variable ``<varName>`` must exist before calling this function.  To set
# it empty initially use `global_null_set()`_.
#
macro(prepend_global_set VARNAME)
  assert_defined(${VARNAME})
  global_set(${VARNAME} ${ARGN} ${${VARNAME}})
endmacro()
