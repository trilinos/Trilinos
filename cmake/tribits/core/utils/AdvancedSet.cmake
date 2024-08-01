# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: advanced_set()
#
# Macro that sets a variable and marks it as advanced (removes boiler-plate
# and duplication).
#
# Usage::
#
#   advanced_set(<varName> [other arguments])
#
# This just calls the built-in commands::
#
#   set(<varName> [other arguments])
#   mark_as_advanced(<varName>)
#
macro(advanced_set VARNAME)
  set(${VARNAME} ${ARGN})
  mark_as_advanced(${VARNAME})
endmacro()
