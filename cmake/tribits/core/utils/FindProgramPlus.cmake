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
