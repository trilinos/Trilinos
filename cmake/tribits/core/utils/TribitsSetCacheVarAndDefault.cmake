# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()


# @MACRO: tribits_advanced_set_cache_var_and_default()
#
# Set an advanced cache variable with a default value (passing in a default
# default value).
#
# Usage::
#
#   tribits_advanced_set_cache_var_and_default(<cacheVarName>  <cacheVarType>
#     <defaultDefaultVal>  <docString>)
#
# If the variable ``<cacheVarName>_DEFAULT`` already exists with a value, that
# is used as the default cache variable.  Otherwise,
# ``<cacheVarName>_DEFAULT`` is set set to ``<defaultDefaultVal>`` first.
#
macro(tribits_advanced_set_cache_var_and_default  cacheVarName  cacheVarType
    defaultDefaultVal  docString
  )
  tribits_set_cache_var_and_default("${cacheVarName}" "${cacheVarType}"
    "${defaultDefaultVal}" "${docString}")
  mark_as_advanced(${cacheVarName})
endmacro()


# @MACRO: tribits_set_cache_var_and_default()
#
# Set a cache variable with a default value (passing in a default default
# value).
#
# Usage::
#
#   tribits_set_cache_var_and_default(<cacheVarName>  <cacheVarType>
#     <defaultDefaultVal>  <docString>)
#
# If the variable ``<cacheVarName>_DEFAULT`` already exists with a value, that
# is used as the default cache variable.  Otherwise,
# ``<cacheVarName>_DEFAULT`` is set set to ``<defaultDefaultVal>`` first.
#
macro(tribits_set_cache_var_and_default  cacheVarName  cacheVarType
    defaultDefaultVal  docString
  )
  if ("${${cacheVarName}_DEFAULT}" STREQUAL "")
    set(${cacheVarName}_DEFAULT "${defaultDefaultVal}")
  endif()
  set(${cacheVarName} "${${cacheVarName}_DEFAULT}"
    CACHE ${cacheVarType}
    "${docString}" )
endmacro()
