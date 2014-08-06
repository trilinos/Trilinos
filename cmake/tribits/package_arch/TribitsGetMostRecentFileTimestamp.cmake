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

INCLUDE(ParseVariableArguments)

#
# @FUNCTION: TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP()
#
# Find the most modified file in a set of base directories and return its
# timestamp.
#
# Usage::
#
#   TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP(
#     BASE_DIRS <dir0> <dir1> ... 
#     [EXCLUDE_REGEXS "<re0>" "<re1>" ... 
#     [MOST_RECENT_TIMESTAMP_OUT  <mostRecentTimestamp>]
#     [MOST_RECENT_FILEPATH_BASE_DIR_OUT <mostRecentFilepathBaseDir>]
#     [MOST_RECENT_RELATIVE_FILEPATH_OUT <mostRecentRelativeFilePath>]
#     [SHOW_MOST_RECENT_FILES]
#     [SHOW_OVERALL_MOST_RECENT_FILE]
#     )
#
FUNCTION(TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP)

    IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DUMP)
      MESSAGE("\nSearching for most modified files in base dirs:")
    ENDIF()
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "BASE_DIRS;EXCLUDE_REGEXS;MOST_RECENT_TIMESTAMP_OUT"
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     ${ARGN}
     )

  IF (PARSE_SHOW_MOST_RECENT_FILES)
    SET(PARSE_SHOW_OVERALL_MOST_RECENT_FILE ON)
  ENDIF()

  #
  # B) Loop over each directory and find the most modified file
  #

  SET(OVERALL_MOST_RECENT_TIMESTAMP "0000000000.0000000000") 
  SET(OVERALL_MOST_RECENT_FILEPATH "") 
  SET(OVERALL_MOST_RECENT_FILEPATH_DIR "") 
  SET(OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE "") 

  FOREACH(BASE_DIR ${PARSE_BASE_DIRS})
    IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DUMP)
      MESSAGE("\nSearching '${BASE_DIR}' ...")
    ENDIF()

    # Get the time stamp and the file name of the most recently modified file
    # in currnet directory.
    EXECUTE_PROCESS(
      WORKING_DIRECTORY "${BASE_DIR}"
      COMMAND find . -type f -printf "'%T@ %p\n'"
      COMMAND sort -n
      COMMAND tail -1
      OUTPUT_VARIABLE MOST_RECENT_TIMESTAMP_AND_FILE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
     # Here, this will return a string with the date and the file name of the
     # form:
     #
     #     ''1407353359.5651538200 ./<relative-dir>/<some-file-name>
     #
     # Note the two chars '' at the beginning.

    IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DUMP)
      PRINT_VAR(MOST_RECENT_TIMESTAMP_AND_FILE)
    ENDIF()

    SPLIT("${MOST_RECENT_TIMESTAMP_AND_FILE}" " "
      MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT)

    # Get the time stamp part and remove the initial "''" chars
    LIST(GET MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT 0
      CURRENT_TIMESTAMP)
    STRING(SUBSTRING "${CURRENT_TIMESTAMP}" 2 -1 CURRENT_TIMESTAMP)

    # Get the relative file path
    LIST(GET MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT 1
      CURRENT_FILEPATH)

    IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DUMP)
      PRINT_VAR(CURRENT_TIMESTAMP)
      PRINT_VAR(CURRENT_FILEPATH)
    ENDIF()
  
    IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
      EXECUTE_PROCESS(
        WORKING_DIRECTORY "${BASE_DIR}"
        COMMAND ls --full-time "${CURRENT_FILEPATH}"
        OUTPUT_VARIABLE HUMAN_READABLE_FILE_AND_TIMESTAMP
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
      IF (PARSE_SHOW_MOST_RECENT_FILES)
        MESSAGE("-- " "Most recent file in '${BASE_DIR}'\n"
          "       ${HUMAN_READABLE_FILE_AND_TIMESTAMP}")
      ENDIF()
    ENDIF()

    IF ("${CURRENT_TIMESTAMP}" GREATER "${OVERALL_MOST_RECENT_TIMESTAMP}")
      SET(OVERALL_MOST_RECENT_TIMESTAMP "${CURRENT_TIMESTAMP}")
      SET(OVERALL_MOST_RECENT_FILEPATH "${BASE_DIR}/${CURRENT_FILEPATH}")
      IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DUMP)
        MESSAGE("    New most recent file path!")
      ENDIF()
      IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
        SET(OVERALL_MOST_RECENT_FILEPATH_DIR "${BASE_DIR}") 
        SET(OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE
          "${HUMAN_READABLE_FILE_AND_TIMESTAMP}") 
      ENDIF()
    ENDIF()

  ENDFOREACH()

  IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DUMP)
    PRINT_VAR(OVERALL_MOST_RECENT_TIMESTAMP)
    PRINT_VAR(OVERALL_MOST_RECENT_FILEPATH)
  ENDIF()

  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    MESSAGE("-- " "Overall most recent file is in:\n"
      "      ${OVERALL_MOST_RECENT_FILEPATH_DIR}\n"
      "       ${HUMAN_READABLE_FILE_AND_TIMESTAMP}")
  ENDIF()

  SET(${PARSE_MOST_RECENT_TIMESTAMP_OUT} ${OVERALL_MOST_RECENT_TIMESTAMP}
    PARENT_SCOPE)

ENDFUNCTION()


#
# @FUNCTION: TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP()
#
# Find the most modified file in a set of base directories and return its
# timestamp.
#
# Usage::
#
#   TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP(
#     SOURCE_BASE_DIRS <dir0> <dir1> ... 
#     [MOST_RECENT_TIMESTAMP_OUT  <mostRecentTimestamp>]
#     [MOST_RECENT_FILEPATH_BASE_DIR_OUT <mostRecentFilepathBaseDir>]
#     [MOST_RECENT_RELATIVE_FILEPATH_OUT <mostRecentRelativeFilePath>]
#     [SHOW_MOST_RECENT_FILES]
#     [SHOW_OVERALL_MOST_RECENT_FILE]
#     )
#
FUNCTION(TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP)
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "SOURCE_BASE_DIRS;MOST_RECENT_TIMESTAMP_OUT"
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     ${ARGN}
     )

  #
  # B) Call the function TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP()
  #

  SET(ARGS)
  APPEND_SET(ARGS BASE_DIRS ${PARSE_SOURCE_BASE_DIRS})
  APPEND_SET(ARGS MOST_RECENT_TIMESTAMP_OUT MOST_RECENT_TIMESTAMP)
  IF (PARSE_SHOW_MOST_RECENT_FILES)
    APPEND_SET(ARGS SHOW_MOST_RECENT_FILES)
  ENDIF()
  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    APPEND_SET(ARGS SHOW_OVERALL_MOST_RECENT_FILES)
  ENDIF()

  #PRINT_VAR(ARGS)
  TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP(${ARGS})
  #PRINT_VAR(MOST_RECENT_TIMESTAMP)

  SET(${PARSE_MOST_RECENT_TIMESTAMP_OUT} ${MOST_RECENT_TIMESTAMP}
    PARENT_SCOPE)

ENDFUNCTION()
