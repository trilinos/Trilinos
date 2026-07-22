# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsDeprecatedHelpers)


macro(parse_arguments_deprecated_warning)
  tribits_deprecated_command(parse_arguments
    MESSAGE "Use cmake_parse_arguments() instead.")
endmacro()

parse_arguments_deprecated_warning()

# Set PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED to TRUE to see output from parsing.

function(parse_arguments_dump_output  OUTPUT_STR)
  if (PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED)
    message("${OUTPUT_STR}")
  endif()
endfunction()


# @MACRO: parse_arguments()
#
# Parse a set of macro/function input arguments into different lists.  This
# allows the easy implementation of keyword-based user-defined macros and
# functions.
#
# Usage::
#
#   parse_arguments(
#     <prefix>  <argNamesList>  <optionNamesList>
#     <inputArgsList>
#     )
#
# Arguments to this macro are:
#
#   ``<prefix>``
#
#     Prefix ``<prefix>_`` added the list and option variables created listed
#     in ``<argNamesList>`` and ``<optionNamesList>``.
#
#   ``<argNamesList>``
#
#     Quoted array of list arguments (e.g. ``"<argName0>;<argName1>;..."``).
#     For each variable name ``<argNamei>``, a local variable will be created
#     in the current scope with the name ``<prefix>_<argNamei>`` which gives
#     the list of variables parsed out of ``<inputArgsList>``.
#
#   ``<optionNamesList>``
#
#     Quoted array of list options (e.g. ``"<optName0>;<optName1>;..."``)
#     typically pass in as ``${ARGN}`` in the outer function/macro.  For each
#     variable name ``<optNamei>``, a local variable will be created in the
#     current scope with the name ``<prefix>_<optNamei>`` that is either set
#     to ``TRUE`` or ``FALSE`` depending if ``<optNamei>`` appears in
#     ``<inputArgsList>`` or not.
#
#   ``<inputArgsList>``
#
#     List of arguments keyword-based arguments passed in through the outer
#     macro or function to be parsed out into the different argument and
#     option lists.
#
# What this macro does is very simple yet very useful.  What it does is to
# allow one to create one's own user-defined keyword-based macros and
# functions like is used by some built-in CMake commands.
#
# For example, consider the following user-defined macro that uses both
# positional and keyword-based arguments using ``parse_arguments()``::
#
#   macro(parse_special_vars  BASE_NAME)
#
#     parse_arguments(
#       #prefix
#       ${BASE_NAME}
#       #lists
#       "ARG0;ARG1;ARG2"
#       #options
#       "OPT0;OPT1"
#       ${ARGN}
#       )
#
#   endmacro()
#
# Calling this macro as::
#
#   parse_special_vars(MyVar ARG0 a b ARG2 c OPT1)
#
# sets the following variables in the current scope::
#
#   MyVar_ARG0="a;b"
#   MyVar_ARG1=""
#   MyVar_ARG2="c"
#   MyVar_OPT0="FALSE"
#   MyVar_OPT1="TRUE"
#
# This allows one to define user-defined macros and functions that have a
# mixture of positional arguments and keyword-based arguments like one can do
# in other languages.  The keyword-based arguments can be passed in any order
# and those that are missing are empty (or false for options) by default.
#
# Any initial arguments that are not recognized as ``<argNamesList>`` or
# ``<optionNamesList>`` keyword arguments will be put into the local variable
# ``<prefix>_DEFAULT_ARGS``.  If no arguments in ``<inputArgsList>``
# (typically ``${ARGN}``) match any in ``<argNamesList>``, then all non-option
# arguments are put into ``<prefix>_DEFAULT_ARGS``.  For example, if one
# passes in::
#
#   parse_special_vars(MyVar ARG5 a b c)
#
# you will get::
#
#   MyVar_DEFAULT_ARGS="ARG5;a;b;c"
#   MyVar_ARG0=""
#   MyVar_ARG1=""
#   MyVar_ARG2=""
#   MyVar_OPT0="FALSE"
#   MyVar_OPT1="FALSE"
#
# Multiple occurrences of keyword arguments in ``<inputArgsList>`` is allowed
# but only the last one listed will be recorded.  For example, if one calls::
#
#   parse_special_vars(MyVar ARG1 a b ARG1 c)
#
# then this will set::
#
#   MyVar_ARG0=""
#   MyVar_ARG1="c"
#   MyVar_ARG2=""
#   MyVar_OPT0="FALSE"
#   MyVar_OPT1="FALSE"
#
# This is actually consistent with the way that most argument list parsers
# behave with respect to multiple instances of the same argument so hopefully
# this will not be a surprise to anyone.
#
# If one puts an option keyword in the middle of a keyword argument list, the
# option keyword will get pulled out of the list.  For example, if one calls::
#
#   parse_special_vars(MyVar ARG0 a OPT0 c)
#
# then this will set::
#
#   MyVar_ARG0="a;c"
#   MyVar_ARG1=""
#   MyVar_ARG2=""
#   MyVar_OPT0="TRUE"
#   MyVar_OPT1="FALSE"
#
# This is confusing behavior so users would be smart not to mix option
# arguments inside of list arguments.
#
# If ``PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED`` is set to ``TRUE``, then a bunch
# of detailed debug info will be printed.  This should only be used in the
# most desperate of debug situations because it will print a *lot* of output!
#
# **PERFORMANCE:** This function will scale as::
#
#   o( (len(<argNamesList>) * len(<optionNamesList>)) * len(<inputArgsList>) )
#
# Therefore, this could scale very badly for large sets of argument and option
# names and input argument list names.
#
macro(parse_arguments prefix arg_names option_names)
   
  parse_arguments_deprecated_warning()
 
  parse_arguments_dump_output("PARSE_ARGUMENTS: prefix='${prefix}'")
  parse_arguments_dump_output("PARSE_ARGUMENTS: arg_names='${arg_names}'")
  parse_arguments_dump_output("PARSE_ARGUMENTS: option_names='${option_names}'")
  parse_arguments_dump_output("PARSE_ARGUMENTS: ARGN='${ARGN}'")

  foreach(arg_name ${arg_names})
    set(${prefix}_${arg_name})
  endforeach()

  foreach(option ${option_names})
    set(${prefix}_${option} FALSE)
  endforeach()

  set(DEFAULT_ARGS)
  set(current_arg_name DEFAULT_ARGS)
  set(current_arg_list)

  foreach(arg ${ARGN})
    set(larg_names ${arg_names})
    list(FIND larg_names "${arg}" is_arg_name)
    if (is_arg_name GREATER -1)
      set(${prefix}_${current_arg_name} "${current_arg_list}")
      parse_arguments_dump_output("PARSE_ARGUMENTS: ${prefix}_${current_arg_name} = '${${prefix}_${current_arg_name}}'" )
      set(current_arg_name "${arg}")
      set(current_arg_list)
    else()
      set(loption_names "${option_names}")
      list(FIND loption_names "${arg}" is_option)
      if (is_option GREATER -1)
        set(${prefix}_${arg} TRUE)
        parse_arguments_dump_output( "PARSE_ARGUMENTS: ${prefix}_${arg} = '${${prefix}_${arg}}'" )
      else()
        list(APPEND current_arg_list "${arg}")
      endif()
    endif()
  endforeach()

  set(${prefix}_${current_arg_name} "${current_arg_list}")
  parse_arguments_dump_output( "PARSE_ARGUMENTS: ${prefix}_${current_arg_name} = '${${prefix}_${current_arg_name}}'" )

endmacro()

# NOTE: If the above function turns out to be a performance bottle neck, there
# are a few things that could be done to improve performance.  One thing you
# could do is replace the o(len(arg_names)) and o(len(option_names)) lookup
# with o(1) lookups by creating CMake variables of the name
# ${OUTER_FUNC_NAME}_arg_<argNamei> and then just look of that variable exists
# or not.  That should use a hash function.  That might actually slow things
# down for short lists however so we would have to measure, measure,
# measure. I would have to pass in the function/macro name to disambiguate
# the variable names.  It would really be better if CMake would provide a
# sorted list find operation.  That would make this much faster for large
# numbers of argument and option names.
