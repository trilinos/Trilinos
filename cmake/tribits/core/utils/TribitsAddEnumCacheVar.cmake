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

include(TribitsParseArgumentsHelpers)


# @FUNCTION: tribits_add_enum_cache_var()
#
# Set up a string cache variable that must match a fixed set of values
# (i.e. an enum) and assert that it matches those values.
#
# Usage::
#
#   tribits_add_enum_cache_var(<cacheVarName>
#     DEFAULT_VAL <defaultVal>
#     DOC_STRING "<docString>"
#     ALLOWED_STRINGS_LIST "<val0>" "<val1>" ...
#     [IS_ADVANCED]
#     )
#
# On output, ``<cacheVarName>`` will be set to the list of paths 
#
function(tribits_add_enum_cache_var  cacheVarName  defaultVal  docString
    isAdvanced  allowedStringsListName
  )
  # Parse input arguments
  set(argOneValArgKeywords  DEFAULT_VAL  DOC_STRING)
  set(argMultiValArgKeywords  ALLOWED_STRINGS_LIST)
  cmake_parse_arguments(PARSE_ARGV  1  PREFIX
    "IS_ADVANCED"   # options
    ""   # one_value_keywords
    "${argOneValArgKeywords};${argMultiValArgKeywords}"   # multi_value_keywords
    )
  tribits_check_for_unparsed_arguments(PREFIX)
  tribits_assert_parse_arg_one_value(PREFIX ${argOneValArgKeywords}) 
  tribits_assert_parse_arg_one_or_more_values(PREFIX ${argMultiValArgKeywords}) 
  # Human readable list of allowed values: '<val0>', '<val1>', ...
  string(REPLACE ";" "', '" validStringValuesListStr "'${PREFIX_ALLOWED_STRINGS_LIST}'")
  # Set cache var
  set(${cacheVarName}  ${PREFIX_DEFAULT_VAL}  CACHE  STRING
    "${PREFIX_DOC_STRING}.  Valid values: ${validStringValuesListStr} (default '${PREFIX_DEFAULT_VAL}')")
  if (PREFIX_IS_ADVANCED)
    mark_as_advanced(${cacheVarName})
  endif()
  set_property(CACHE  ${cacheVarName}  PROPERTY  STRINGS
    ${PREFIX_ALLOWED_STRINGS_LIST} )
  # Assert in list of allowed strings
  if (NOT  ${cacheVarName}  IN_LIST  PREFIX_ALLOWED_STRINGS_LIST)
    message(FATAL_ERROR "Error, the cache var ${cacheVarName} with value"
      " '${${cacheVarName}}' is not in the list of allowed values:"
      " ${validStringValuesListStr} (default '${PREFIX_DEFAULT_VAL}')"  )
  endif()
endfunction()
