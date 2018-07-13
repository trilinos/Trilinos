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

INCLUDE(CMakeParseArguments)
INCLUDE(MultilineSet)
INCLUDE(ConcatStrings)


#
# @FUNCTION: COMBINED_OPTION()
#
# Set up a ``BOOL`` cache variable (i.e. an option) based on a set of
# dependent options.
#
# Usage::
#
#   COMBINED_OPTION( <combinedOptionName>
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
FUNCTION(COMBINED_OPTION  COMBINED_OPTION_NAME)

  CMAKE_PARSE_ARGUMENTS(
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

  TRIBITS_CHECK_FOR_UNPARSED_ARGUMENTS()

  # ToDo: Assert that the right input was passed in!

  SET(DEFAULT_VAL ON)
  FOREACH( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES} )
    IF (NOT ${DEP_OPTION_NAME})
      SET(DEFAULT_VAL OFF)
    ENDIF()
  ENDFOREACH()

  CONCAT_STRINGS(DOCSTR ${PARSE_DOCSTR})

  OPTION(${COMBINED_OPTION_NAME} "${DOCSTR}"
    ${DEFAULT_VAL} )

  #PRINT_VAR(${COMBINED_OPTION_NAME})

  # Determine if the combined option was turned on by the individual options
  # are not turned on as well.
  SET(ALL_ENABLED TRUE)
  IF (${COMBINED_OPTION_NAME})
    FOREACH( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES})
      #PRINT_VAR(${DEP_OPTION_NAME})
      IF (NOT ${DEP_OPTION_NAME})
        SET(ALL_ENABLED FALSE)
        BREAK()
      ENDIF()
    ENDFOREACH()
  ENDIF()

  # Print out detailed error message if the combined option was enabled but
  # the dependent options were not.
  IF (NOT ALL_ENABLED)

    SET(OPTION_NAMES "")
    SET(OPTION_NAMES_AND_VALUES "")
    FOREACH( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES})
      IF (NOT OPTION_NAMES)
        SET(OPTION_NAMES "${DEP_OPTION_NAME}")
      ELSE()
        SET(OPTION_NAMES "${OPTION_NAMES}, ${DEP_OPTION_NAME}")
      ENDIF()
      SET(OPTION_NAMES_AND_VALUES
        "${OPTION_NAMES_AND_VALUES}  ${DEP_OPTION_NAME}='${${DEP_OPTION_NAME}}'\n")
    ENDFOREACH()

    MESSAGE(FATAL_ERROR
      "Error: you can not enable the option ${COMBINED_OPTION_NAME} unless"
      " you also enable the options ${OPTION_NAMES}.  The current option"
      " values are:\n${OPTION_NAMES_AND_VALUES}" )

  ENDIF()

ENDFUNCTION()
