# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: advanced_option()
#
# Macro that sets an option and marks it as advanced (removes boiler-plate and
# duplication).
#
# Usage::
#
#   advanced_option(<varName> [other arguments])
#
# This just calls the built-in CMake commands::
#
#   option(<varName> [other arguments])
#   mark_as_advanced(<varName>)
#
macro(advanced_option VARNAME)
  option(${VARNAME} ${ARGN})
  mark_as_advanced(${VARNAME})
endmacro()
