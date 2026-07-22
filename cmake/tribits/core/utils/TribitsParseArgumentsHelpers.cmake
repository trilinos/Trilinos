# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


################################################################################
#
# This module contains function to aid with the usage of
# cmake_parse_arguments().
#
################################################################################


include("${CMAKE_CURRENT_LIST_DIR}/MessageWrapper.cmake")


# @FUNCTION: tribits_check_for_unparsed_arguments()
#
# Check to see if there are unparsed arguments after calling
# ``cmake_parse_arguments()`` or ``tribits_parse_arguments_from_list()``
#
# Usage::
#
#   tribits_check_for_unparsed_arguments([<prefix>])
#
# If `<prefix>` is not given, it is assumed to be `PARSE`.
#
function(tribits_check_for_unparsed_arguments)

  if ("${ARGC}" GREATER 1)
    message(FATAL_ERROR
      "ERROR tribits_check_for_unparsed_arguments() passed arguments '${ARGV}' but only accepts 0 or 1 arguments (for <prefix>)")
  endif()

  set(prefix "PARSE")
  foreach(arg ${ARGV})
    set(prefix "${arg}")
  endforeach()

  if( NOT "${${prefix}_UNPARSED_ARGUMENTS}" STREQUAL "")
    message_wrapper(
      ${${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS}
      "Arguments passed in unrecognized.  ${prefix}_UNPARSED_ARGUMENTS = '${${prefix}_UNPARSED_ARGUMENTS}'"
      )
  endif()

endfunction()


# @MACRO: tribits_assert_parse_arg_one_or_more_values()
#
# Assert that a set of parse arguments have at least one value
#
# Usage::
#
#   tribits_assert_parse_arg_one_or_more_values(<prefix> <argname0> <argname1> ...)
#
macro(tribits_assert_parse_arg_one_or_more_values  PREFIX)
  foreach(ARGNAME ${ARGN})
    set(PREFIX_ARGNAME "${PREFIX}_${ARGNAME}")
    list( LENGTH ${PREFIX_ARGNAME} ARG_NUM_VALS )
    if (ARG_NUM_VALS LESS 1)
      message_wrapper(FATAL_ERROR
        "ERROR: ${ARGNAME} must have at least one value!" )
      return()
      # NOTE: The return() is needed in unit testing mode
    endif()
  endforeach()
endmacro()


# @MACRO: tribits_assert_parse_arg_zero_or_one_value()
#
# Assert a set of parse arguments have zero or one value
#
# Usage::
#
#   tribits_assert_parse_arg_zero_or_one_value(<prefix> <argname0> <argname1> ...)
#
macro(tribits_assert_parse_arg_zero_or_one_value  PREFIX)
  foreach(ARGNAME ${ARGN})
    set(PREFIX_ARGNAME "${PREFIX}_${ARGNAME}")
    if (NOT "${${PREFIX_ARGNAME}}" STREQUAL "")
      list( LENGTH ${PREFIX_ARGNAME} ARG_NUM_VALS )
      if (ARG_NUM_VALS GREATER 1)
        message_wrapper(FATAL_ERROR
          "ERROR: ${ARGNAME}='${${PREFIX_ARGNAME}}' can not have more than one value!" )
        return()
        # NOTE: The macro return() is needed in unit testing mode
      endif()
    endif()
  endforeach()
endmacro()


# @MACRO: tribits_assert_parse_arg_one_value()
#
# Assert that a set of parse arguments have exactly one value
#
# Usage::
#
#   tribits_assert_parse_arg_one_value(<prefix> <argname0> <argname1> ...)
#
macro(tribits_assert_parse_arg_one_value  PREFIX)
  foreach(ARGNAME ${ARGN})
    set(PREFIX_ARGNAME "${PREFIX}_${ARGNAME}")
    if (NOT "${${PREFIX_ARGNAME}}" STREQUAL "")
      list( LENGTH ${PREFIX_ARGNAME} ARG_NUM_VALS )
      if (NOT ARG_NUM_VALS EQUAL 1)
        message_wrapper(FATAL_ERROR
          "ERROR: ${ARGNAME}='${${PREFIX_ARGNAME}}' Must have exactly one value!" )
        return()
        # NOTE: The macro return() is needed in unit testing mode
      endif()
    endif()
  endforeach()
endmacro()


# NOTE: Above, we use macros for the assert functions with returns in unit
# test mode so that it will abort the calling function these are called from!
