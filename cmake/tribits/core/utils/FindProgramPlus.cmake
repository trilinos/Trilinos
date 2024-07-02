# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(CMakeParseArguments)
include(PrintNonemptyVar)

#
# Finds the absolute path for a program given optionally just the program name
#
#  Input keyword arguments:
#
#    NAMES name1 name2 ...
#     List of default names that will be search
#
#    PATHS path1 path2 ...
#      List of paths
#
#    DOC docStr
#      Documentation string
#

function(find_program_plus PROG_VAR)

  cmake_parse_arguments(
    #prefix
    PARSE
    #options
    ""
    #one_value_keywords
    ""
    #multi_value_keywords
    "NAMES;PATHS;DOC"
    ${ARGN}
    )

  tribits_check_for_unparsed_arguments()

  print_nonempty_var(${PROG_VAR})

  if (IS_ABSOLUTE ${PROG_VAR})
    #message(STATUS "Is Absolute")
    set(NAMES_ARGS ${PARSE_NAMES})
  else()
    #message(STATUS "Is Not Absolute")
    set(NAMES_ARGS ${${PROG_VAR}} ${PARSE_NAMES})
    set(${PROG_VAR} "${PROG_VAR}-NOTFOUND" CACHE FILEPATH "" FORCE)
  endif()
  #print_var(NAMES_ARGS)

  set(DOC "${PARSE_DOC}  Can be full path or just exec name.")

  # Look for program in given paths first!
  find_program( ${PROG_VAR}
    NAMES ${NAMES_ARGS}
    PATHS ${PARSE_PATHS}
    DOC ${DOC}
    NO_DEFAULT_PATH
    )
  find_program( ${PROG_VAR}
    NAMES ${NAMES_ARGS}
    DOC ${DOC}
    )
  mark_as_advanced(${PROG_VAR})

  print_nonempty_var(${PROG_VAR})

endfunction()
