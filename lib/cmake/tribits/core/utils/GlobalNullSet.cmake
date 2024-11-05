# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: global_null_set()
#
# Set a variable as a null internal global (cache) variable (removes
# boiler-plate code).
#
# Usage::
#
#   global_null_set(<varName>)
#
# This just calls::
#
#   set(<varName> "" CACHE INTERNAL "")
#
# This avoid problems with misspelling ``CACHE``.
#
macro(global_null_set VARNAME)
  set(${VARNAME} "" CACHE INTERNAL "")
endmacro()
