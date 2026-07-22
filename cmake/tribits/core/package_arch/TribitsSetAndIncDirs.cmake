# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: tribits_set_and_inc_dirs()
#
# Set a variable to an include directory and call
# `tribits_include_directories()`_ (removes boiler-plate code).
#
# Usage:
#
#   tribits_set_and_inc_dirs(<dirVarName> <includeDir>)
#
# On output, this sets ``<dirVarName>`` to ``<includeDir>`` in the local scope
# and calls ``tribits_include_directories(<includeDir>)``.
#
macro(tribits_set_and_inc_dirs  dirVarName  includeDir)
  set(${dirVarName} ${includeDir})
  tribits_include_directories(${${dirVarName}})
endmacro()


# Deprecated!  Use tribits_set_and_inc_dirs() instead!
#
macro(set_and_inc_dirs  DIR_VAR_NAME  INCLUDE_DIR)
  tribits_deprecated_command(set_and_inc_dirs
    MESSAGE "Use tribits_set_and_inc_dirs() instead." )
  set(${DIR_VAR_NAME} ${INCLUDE_DIR})
  tribits_include_directories(${${DIR_VAR_NAME}})
endmacro()
