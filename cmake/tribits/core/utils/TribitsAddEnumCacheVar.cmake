# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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
