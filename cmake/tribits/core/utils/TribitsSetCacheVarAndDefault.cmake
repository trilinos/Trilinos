# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
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
