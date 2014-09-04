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


# Set PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED to TRUE to see output from parsing.

FUNCTION(PARSE_ARGUMENTS_DUMP_OUTPUT  OUTPUT_STR)
  IF (PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED)
    MESSAGE("${OUTPUT_STR}")
  ENDIF()
ENDFUNCTION()

#
# @MACRO: PARSE_ARGUMENTS()
#
# Parse a set of macro/function input arguments into different lists.  This
# allows the easy implementation of keyword-based user-defined macros and
# functions.
#
# Usage::
#
#   PARSE_ARGUMENTS(
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
# positional and keyword-based arguments using ``PARSE_ARGUMENTS()``::
#
#   MACRO(PARSE_SPECIAL_VARS  BASE_NAME)
#
#     PARSE_ARGUMENTS(
#       #prefix
#       ${BASE_NAME}
#       #lists
#       "ARG0;ARG1;ARG2"
#       #options
#       "OPT0;OPT1"
#       ${ARGN}
#       )
#
#   ENDMACRO()
#
# Calling this macro as::
#
#   PARSE_SPECIAL_VARS(MyVar ARG0 a b ARG2 c OPT1)
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
#   PARSE_SPECIAL_VARS(MyVar ARG5 a b c)
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
#   PARSE_SPECIAL_VARS(MyVar ARG1 a b ARG1 c)
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
#   PARSE_SPECIAL_VARS(MyVar ARG0 a OPT0 c)
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
#   O( (len(<argNamesList>) * len(<optionNamesList>)) * len(<inputArgsList>) )
#
# Therefore, this could scale very badly for large sets of argument and option
# names and input argument list names.
#
MACRO(PARSE_ARGUMENTS prefix arg_names option_names)

  PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: prefix='${prefix}'")
  PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: arg_names='${arg_names}'")
  PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: option_names='${option_names}'")
  PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: ARGN='${ARGN}'")

  FOREACH(arg_name ${arg_names})    
    SET(${prefix}_${arg_name})
  ENDFOREACH()

  FOREACH(option ${option_names})
    SET(${prefix}_${option} FALSE)
  ENDFOREACH()

  SET(DEFAULT_ARGS)
  SET(current_arg_name DEFAULT_ARGS)
  SET(current_arg_list)

  FOREACH(arg ${ARGN})            
    SET(larg_names ${arg_names})    
    LIST(FIND larg_names "${arg}" is_arg_name)                   
    IF (is_arg_name GREATER -1)
      SET(${prefix}_${current_arg_name} "${current_arg_list}")
      PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: ${prefix}_${current_arg_name} = '${${prefix}_${current_arg_name}}'" )
      SET(current_arg_name "${arg}")
      SET(current_arg_list)
    ELSE()
      SET(loption_names "${option_names}")
      LIST(FIND loption_names "${arg}" is_option)            
      IF (is_option GREATER -1)
        SET(${prefix}_${arg} TRUE)
        PARSE_ARGUMENTS_DUMP_OUTPUT( "PARSE_ARGUMENTS: ${prefix}_${arg} = '${${prefix}_${arg}}'" )
      ELSE()
        LIST(APPEND current_arg_list "${arg}")
      ENDIF()
    ENDIF()
  ENDFOREACH()

  SET(${prefix}_${current_arg_name} "${current_arg_list}")
  PARSE_ARGUMENTS_DUMP_OUTPUT( "PARSE_ARGUMENTS: ${prefix}_${current_arg_name} = '${${prefix}_${current_arg_name}}'" )

ENDMACRO()

# NOTE: If the above function turns out to be a performance bottle neck, there
# are a few things that could be done to improve performance.  One thing you
# could do is repalce the O(len(arg_names)) and O(len(option_names)) lookups
# with O(1) lookups by creating CMake varibles of the name
# ${OUTER_FUNC_NAME}_arg_<argNamei> and then just look of that varible exists
# or not.  That should use a hash function.  That might actually slow things
# down for short lists however so we would have to measure, measure,
# measure. I we would have to pass in the function/macro name to disabiguate
# the varible names.  It would really be better if cmake would provide a
# sorted list find operation.  That would make this much faster for large
# numbers of argument and option names.
