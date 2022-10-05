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
include(MultilineSet)
include(ConcatStrings)


# @FUNCTION: combined_option()
#
# Set up a ``BOOL`` cache variable (i.e. an option) based on a set of
# dependent options.
#
# Usage::
#
#   combined_option( <combinedOptionName>
#     DEP_OPTIONS_NAMES <depOpName0> <depOptName1> ...
#     DOCSTR "<docstr0>" "<docstr1>" ...
#     )
#
# This sets up a ``BOOL`` cache variable ``<combinedOptionName>`` which is
# defaulted to ``ON`` if all of the listed dependent option variables
# ``<depOpName0>``, ``<depOptName1>``, ... are all ``ON``.  However, if
# ``<combinedOptionName>`` is set to ``ON`` by the user and not all of the
# dependent option variables are also ``ON``, then this results in a fatal
# error and all processing stops.
#
# This is used by a CMake project to automatically turn on a feature that
# requires a set of other features (when they are all enabled) but allows a
# user to disable the feature if desired.
#
function(combined_option  COMBINED_OPTION_NAME)

  cmake_parse_arguments(
    #prefix
    PARSE
    #options
    ""
    #one_value_keywords
    ""
    #multi_value_keywords
    "DEP_OPTIONS_NAMES;DOCSTR"
    ${ARGN}
    )

  tribits_check_for_unparsed_arguments()

  # ToDo: Assert that the right input was passed in!

  set(DEFAULT_VAL ON)
  foreach( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES} )
    if (NOT ${DEP_OPTION_NAME})
      set(DEFAULT_VAL OFF)
    endif()
  endforeach()

  concat_strings(DOCSTR ${PARSE_DOCSTR})

  option(${COMBINED_OPTION_NAME} "${DOCSTR}"
    ${DEFAULT_VAL} )

  #print_var(${COMBINED_OPTION_NAME})

  # Determine if the combined option was turned on by the individual options
  # are not turned on as well.
  set(ALL_ENABLED TRUE)
  if (${COMBINED_OPTION_NAME})
    foreach( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES})
      #print_var(${DEP_OPTION_NAME})
      if (NOT ${DEP_OPTION_NAME})
        set(ALL_ENABLED FALSE)
        break()
      endif()
    endforeach()
  endif()

  # Print out detailed error message if the combined option was enabled but
  # the dependent options were not.
  if (NOT ALL_ENABLED)

    set(OPTION_NAMES "")
    set(OPTION_NAMES_AND_VALUES "")
    foreach( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES})
      if (NOT OPTION_NAMES)
        set(OPTION_NAMES "${DEP_OPTION_NAME}")
      else()
        set(OPTION_NAMES "${OPTION_NAMES}, ${DEP_OPTION_NAME}")
      endif()
      set(OPTION_NAMES_AND_VALUES
        "${OPTION_NAMES_AND_VALUES}  ${DEP_OPTION_NAME}='${${DEP_OPTION_NAME}}'\n")
    endforeach()

    message(FATAL_ERROR
      "Error: you can not enable the option ${COMBINED_OPTION_NAME} unless"
      " you also enable the options ${OPTION_NAMES}.  The current option"
      " values are:\n${OPTION_NAMES_AND_VALUES}" )

  endif()

endfunction()
