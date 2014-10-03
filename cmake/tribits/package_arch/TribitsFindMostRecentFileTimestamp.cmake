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
#     [BASE_BASE_DIR <dir>]
#     [EXCLUDE_REGEXES "<re0>" "<re1>" ... 
#     [SHOW_MOST_RECENT_FILES]
#     [SHOW_OVERALL_MOST_RECENT_FILE]
#     [MOST_RECENT_TIMESTAMP_OUT  <mostRecentTimestamp>]
#     [MOST_RECENT_FILEPATH_BASE_DIR_OUT <mostRecentFilepathBaseDir>]
#     [MOST_RECENT_RELATIVE_FILEPATH_OUT <mostRecentRelativeFilePath>]
#     )
#
# **Arguments:**
#
#   ``BASE_DIRS <dir0> <dir1> ...``
#
#     Gives the absolute base directory paths that will be searched for the
#     most recently modified files, as described above.
#
#   ``BASE_BASE_DIR <dir>```
#
#     Absolute path for which to print file paths relative to.  This makes
#     outputting less verbose and easier to read (optional).
#
#   ``EXCLUDE_REGEXES "<re0>" "<re1>" ...``
#
#     Gives the regular expressions that are used to exclude files from
#     consideration.  Each "<rei>" regex is used with a `grep -v "<rei>"`
#     filter to exclude files before sorting by time stamp.
#
#   ``SHOW_MOST_RECENT_FILES``
#
#     If specified, then the most recently modified file for each individual
#     directory ``<dir0>``, ``<dir1``, ... will be printed the STDOUT.
#     Setting this implies ``SHOW_OVERALL_MOST_RECENT_FILE``.
#
#   ``SHOW_OVERALL_MOST_RECENT_FILE``
#
#     If specified, then only the most recent modified file over all of the
#     individual directories is printed to STDOUT.
#
#   ``MOST_RECENT_TIMESTAMP_OUT <mostRecentTimestamp>``
#
#      On output, the variable `<mostRecentTimestamp>` is set that gives the
#      timestamp of the most recently modified file over all the directories.
#      This number is given as the number of seconds since Jan. 1, 1970, 00:00
#      GMT.
#
#   ``MOST_RECENT_FILEPATH_BASE_DIR_OUT <mostRecentFilepathBaseDir>``
#
#     On output, the variable `<mostRecentFilepathBaseDir>` gives absolute base
#     directory of the file with the most recent timestamp over all
#     directories.
#
#   ``MOST_RECENT_RELATIVE_FILEPATH_OUT <mostRecentRelativeFilePath>``
#
#     On output, the variable `<mostRecentFilepathBaseDir>` gives the file
#     name with relative path to the file with the most recent timestamp over
#     all directories.
#
# **Description:**
#
# This function uses the Linux/Unix command::
#     
#     $ find . -type f -printf '%T@ %p\n'
#         | grep -v "<re0>" | grep -v "<re1>" | ... \
#         | sort -n | tail -1
#
# to return the most recent file in each listed directory <dir0>, <dir1>, etc.
# It then determines the most recently modified file over all of the
# directories and prints and returns in the variables `<mostRecentTimestamp>`,
# `<mostRecentFilepathBaseDir>`, and `<mostRecentRelativeFilePath>`.
#
FUNCTION(TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP)

    IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
      MESSAGE("\nSearching for most modified files in base dirs:")
    ENDIF()
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "BASE_DIRS;BASE_BASE_DIR;EXCLUDE_REGEXES;MOST_RECENT_TIMESTAMP_OUT;MOST_RECENT_FILEPATH_BASE_DIR_OUT;MOST_RECENT_RELATIVE_FILEPATH_OUT"
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
  SET(OVERALL_MOST_RECENT_RELATEIVE_FILEPATH_DIR "")
  SET(OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE "") 

  FOREACH(BASE_DIR ${PARSE_BASE_DIRS})

    IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
      MESSAGE("\nSearching '${BASE_DIR}' ...")
    ENDIF()

    IF (IS_DIRECTORY "${BASE_DIR}")

      # Build up commands for grep -v
      SET(GREP_V_COMMANDS)
      FOREACH(EXCLUDE_REGEX ${PARSE_EXCLUDE_REGEXES})
        APPEND_SET(GREP_V_COMMANDS COMMAND grep -v "${EXCLUDE_REGEX}")
      ENDFOREACH()
  
      # Get the time stamp and the file name of the most recently modified file
      # in currnet directory.
      EXECUTE_PROCESS(
        WORKING_DIRECTORY "${BASE_DIR}"
        COMMAND find . -type f -printf "%T@ %p\n"
        ${GREP_V_COMMANDS}
        COMMAND sort -n
        COMMAND tail -1
        OUTPUT_VARIABLE MOST_RECENT_TIMESTAMP_AND_FILE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
       # Here, note that %T@ gives the modification time stamp in seconds since
       # Jan. 1, 1970, 00:00 GMT.  The -printf argument %p gives the file path.
       # This results in the return a string with the modification date (in
       # fractional seconds) and the file name of the form:
       #
       #     1407353359.5651538200 ./<relative-dir>/<some-file-name>
  
      IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
        PRINT_VAR(MOST_RECENT_TIMESTAMP_AND_FILE)
      ENDIF()

      IF (MOST_RECENT_TIMESTAMP_AND_FILE)

        SPLIT("${MOST_RECENT_TIMESTAMP_AND_FILE}" " "
          MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT)
    
        # Get the time stamp part
        LIST(GET MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT 0
          CURRENT_TIMESTAMP)
    
        # Get the relative file path
        LIST(GET MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT 1
          CURRENT_FILEPATH)

        # Get the directory relative to the base base dir
        STRING(REPLACE "${PARSE_BASE_BASE_DIR}/" "./"  RELATIVE_FILEPATH_DIR
          "${BASE_DIR}")

        IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
          PRINT_VAR(CURRENT_TIMESTAMP)
          PRINT_VAR(CURRENT_FILEPATH)
        ENDIF()
      
        IF (PARSE_SHOW_MOST_RECENT_FILES)
          TRIBITS_GET_HUMAN_READABLE_FILE_AND_TIMESTAMP(
            "${BASE_DIR}"  "${CURRENT_FILEPATH}"
            HUMAN_READABLE_FILE_AND_TIMESTAMP
            )
          MESSAGE("-- " "Most recent file in ${RELATIVE_FILEPATH_DIR}/"
            " is ${CURRENT_FILEPATH}\n"
            "    ${HUMAN_READABLE_FILE_AND_TIMESTAMP}")
        ENDIF()
    
        IF ("${CURRENT_TIMESTAMP}" GREATER "${OVERALL_MOST_RECENT_TIMESTAMP}")
          IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
            MESSAGE("    New most recent file path!")
          ENDIF()
          SET(OVERALL_MOST_RECENT_TIMESTAMP "${CURRENT_TIMESTAMP}")
          SET(OVERALL_MOST_RECENT_RELATIVE_FILEPATH "${CURRENT_FILEPATH}")
          SET(OVERALL_MOST_RECENT_FILEPATH_DIR "${BASE_DIR}")
          SET(OVERALL_MOST_RECENT_RELATEIVE_FILEPATH_DIR "${RELATIVE_FILEPATH_DIR}")
          IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
            SET(OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE
              "${HUMAN_READABLE_FILE_AND_TIMESTAMP}") 
          ENDIF()
        ENDIF()

      ENDIF()

    ELSE()

      IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
        MESSAGE("Directory does not exist, skipping ...")
      ENDIF()

    ENDIF()

  ENDFOREACH()

  IF (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
    PRINT_VAR(OVERALL_MOST_RECENT_TIMESTAMP)
    PRINT_VAR(OVERALL_MOST_RECENT_RELATIVE_FILEPATH)
  ENDIF()

  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    IF (OVERALL_MOST_RECENT_FILEPATH_DIR)
      IF (NOT OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE)
        TRIBITS_GET_HUMAN_READABLE_FILE_AND_TIMESTAMP(
          "${OVERALL_MOST_RECENT_FILEPATH_DIR}"
          "${OVERALL_MOST_RECENT_RELATIVE_FILEPATH}"
          OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE
        )
      ENDIF()
      MESSAGE("-- " "Overall most recent modified file is in"
        " ${OVERALL_MOST_RECENT_RELATEIVE_FILEPATH_DIR}/ and is ${OVERALL_MOST_RECENT_RELATIVE_FILEPATH}\n"
        "    ${OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE}")
    ELSE()
      MESSAGE("-- There are no unfiltered files!")
    ENDIF()
  ENDIF()

  SET(${PARSE_MOST_RECENT_TIMESTAMP_OUT}
    ${OVERALL_MOST_RECENT_TIMESTAMP}
    PARENT_SCOPE)

  IF (PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT)
    SET(${PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT}
      ${OVERALL_MOST_RECENT_FILEPATH_DIR}
      PARENT_SCOPE)
  ENDIF()

  IF (PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT)
    SET(${PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT}
      ${OVERALL_MOST_RECENT_RELATIVE_FILEPATH}
      PARENT_SCOPE )
  ENDIF()

ENDFUNCTION()


#
# @FUNCTION: TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP()
#
# Find the most modified source file in a set of base directories and return
# its timestamp.
#
# Usage::
#
#   TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP(
#     SOURCE_BASE_DIRS <dir0> <dir1> ... 
#     [SOURCE_BASE_BASE_DIR <dir>]
#     [SHOW_MOST_RECENT_FILES]
#     [SHOW_OVERALL_MOST_RECENT_FILE]
#     [MOST_RECENT_TIMESTAMP_OUT  <mostRecentTimestamp>]
#     [MOST_RECENT_FILEPATH_BASE_DIR_OUT <mostRecentFilepathBaseDir>]
#     [MOST_RECENT_RELATIVE_FILEPATH_OUT <mostRecentRelativeFilePath>]
#     )
#
# This function just calls `TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP()`_
# passing in a set of basic exclude regexes like ``[.]git/``, ``[.]svn/``,
# etc.  These types of version control files can not possibly directly impact
# the source code.
#
FUNCTION(TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP)
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "SOURCE_BASE_DIRS;SOURCE_BASE_BASE_DIR;MOST_RECENT_TIMESTAMP_OUT;MOST_RECENT_FILEPATH_BASE_DIR_OUT;MOST_RECENT_RELATIVE_FILEPATH_OUT"
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     ${ARGN}
     )

  #
  # B) Call the function TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP()
  #

  SET(FILTER_OUT_SOURCE_FILE_REGEXS
    "/[.]git/"
    )

  SET(VARIABLE_ARGS)
  IF (PARSE_SHOW_MOST_RECENT_FILES)
    APPEND_SET(VARIABLE_ARGS SHOW_MOST_RECENT_FILES)
  ENDIF()
  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    APPEND_SET(VARIABLE_ARGS SHOW_OVERALL_MOST_RECENT_FILES)
  ENDIF()

  #PRINT_VAR(VARIABLE_ARGS)
  TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP(
    BASE_DIRS ${PARSE_SOURCE_BASE_DIRS}
    BASE_BASE_DIR ${PARSE_SOURCE_BASE_BASE_DIR}
    EXCLUDE_REGEXES ${FILTER_OUT_SOURCE_FILE_REGEXS}
    MOST_RECENT_TIMESTAMP_OUT MOST_RECENT_TIMESTAMP
    MOST_RECENT_FILEPATH_BASE_DIR_OUT  MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT  MOST_RECENT_RELATIVE_FILEPATH
    ${VARIABLE_ARGS}
    )
  #PRINT_VAR(MOST_RECENT_TIMESTAMP)

  SET(${PARSE_MOST_RECENT_TIMESTAMP_OUT} ${MOST_RECENT_TIMESTAMP}
    PARENT_SCOPE)

  IF (PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT)
    SET(${PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT} ${MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP}
      PARENT_SCOPE)
  ENDIF()

  IF (PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT)
    SET(${PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT} ${MOST_RECENT_RELATIVE_FILEPATH}
      PARENT_SCOPE )
  ENDIF()

ENDFUNCTION()


#
# @FUNCTION: TRIBITS_FIND_MOST_RECENT_BINARY_FILE_TIMESTAMP()
#
# Find the most modified binary file in a set of base directories and return
# its timestamp.
#
# Usage::
#
#   TRIBITS_FIND_MOST_RECENT_BINARY_FILE_TIMESTAMP(
#     BINARY_BASE_DIRS <dir0> <dir1> ... 
#     [BINARY_BASE_BASE_DIR <dir>]
#     [MOST_RECENT_TIMESTAMP_OUT  <mostRecentTimestamp>]
#     [MOST_RECENT_FILEPATH_BASE_DIR_OUT <mostRecentFilepathBaseDir>]
#     [MOST_RECENT_RELATIVE_FILEPATH_OUT <mostRecentRelativeFilePath>]
#     [SHOW_MOST_RECENT_FILES]
#     [SHOW_OVERALL_MOST_RECENT_FILE]
#     )
#
# This function just calls `TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP()`_
# passing in a set of basic exclude regexes like ``CMakeFiles/``,
# ``[.]cmake$``, and ``/Makefile$``, etc.  These types of files usually don't
# impact the build of downstream software in CMake projects.
#
FUNCTION(TRIBITS_FIND_MOST_RECENT_BINARY_FILE_TIMESTAMP)
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "BINARY_BASE_DIRS;BINARY_BASE_BASE_DIR;MOST_RECENT_TIMESTAMP_OUT;MOST_RECENT_FILEPATH_BASE_DIR_OUT;MOST_RECENT_RELATIVE_FILEPATH_OUT"
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     ${ARGN}
     )

  #
  # B) Define filters for binary files we know are not significant
  #

  SET(FILTER_OUT_BINARY_FILE_REGEXS
    "CMakeFiles/" "[.]cmake$" "/Makefile$"
    )

  #
  # C) Call the function TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP()
  #

  SET(VARIABLE_ARGS)
  IF (PARSE_SHOW_MOST_RECENT_FILES)
    APPEND_SET(VARIABLE_ARGS SHOW_MOST_RECENT_FILES)
  ENDIF()
  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    APPEND_SET(VARIABLE_ARGS SHOW_OVERALL_MOST_RECENT_FILES)
  ENDIF()

  #PRINT_VAR(VARIABLE_ARGS)
  TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP(
    BASE_DIRS ${PARSE_BINARY_BASE_DIRS}
    BASE_BASE_DIR ${PARSE_BINARY_BASE_BASE_DIR}
    EXCLUDE_REGEXES ${FILTER_OUT_BINARY_FILE_REGEXS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT  MOST_RECENT_RELATIVE_FILEPATH
    ${VARIABLE_ARGS}
    )
  #PRINT_VAR(MOST_RECENT_TIMESTAMP)

  SET(${PARSE_MOST_RECENT_TIMESTAMP_OUT} ${MOST_RECENT_TIMESTAMP}
    PARENT_SCOPE)

  IF (PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT)
    SET(${PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT} ${MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP}
      PARENT_SCOPE)
  ENDIF()

  IF (PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT)
    SET(${PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT} ${MOST_RECENT_RELATIVE_FILEPATH}
      PARENT_SCOPE )
  ENDIF()

ENDFUNCTION()


#
# @FUNCTION: TRIBITS_DETERMINE_IF_CURRENT_PACKAGE_NEEDS_REBUILT()
#
# Determine at configure time if any of the upstream dependencies for a
# package require the current package to be rebuilt.
#
# Usage::
#
#   TRIBITS_DETERMINE_IF_CURRENT_PACKAGE_NEEDS_REBUILT(
#     [SHOW_MOST_RECENT_FILES]
#     [SHOW_OVERALL_MOST_RECENT_FILES]
#     CURRENT_PACKAGE_OUT_OF_DATE_OUT <currentPackageOutOfDate>
#     )
#
# **Arguments:**
#
#   ``SHOW_MOST_RECENT_FILES``
#
#     If specified, then the most recently modified file for each individual
#     base source and binary directory searched will be will be printed the
#     STDOUT.  Setting this implies ``SHOW_OVERALL_MOST_RECENT_FILE``.
#
#   ``SHOW_OVERALL_MOST_RECENT_FILE``
#
#     If specified, then only the most recent modified file over all of the
#     individual directories for each category (i.e. one for upstream SE
#     package source dirs, one for upstream SE package binary dirs, one for
#     the package's source dir, and one for the package's own binary dir) is
#     printed to STDOUT.
#
#   ``CURRENT_PACKAGE_OUT_OF_DATE_OUT <currentPackageOutOfDate>``
#
#     On output, the local variable ``<currentPackageOutOfDate>`` will be set
#     to ``TRUE`` if any of the upstream most modified files are more recent
#     than the most modified file in the package's binary directory.
#     Otherwise, this variable is set to ``FALSE``.
#
# **Description:**
#
# This function is designed to help take an externally configured and built
# piece of software (that generates libraries) and wrap it as a TriBITS
# package or subpackage.  This function uses the lower-level functions:
#
# * `TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP()`_ 
# * `TRIBITS_FIND_MOST_RECENT_BINARY_FILE_TIMESTAMP()`_
#
# to determine the most recent modified files in the upstream TriBITS SE
# packages' source and binary directories as well as the most recent source
# file for the current package.  It then compares these timestamps to the most
# recent binary file timestamp in this package's binary directory.  If any of
# these three files are more recent than this package's most recent binary
# file, then the output variable ``<currentPackageOutOfDate>`` is set to
# ``TRUE``.  Otherwise, it is set to ``FALSE``. 
#
# NOTE: The source and binary directories for full packages are searched, not
# individual subpackage dirs.  This is to reduce the number of dirs searched.
# This will, however, result in changes in non-dependent subpackages being
# considered as well.
#
# See the demonstration of the usage of this function in the ``WrapExternal``
# package in `TribitsExampleProject`_.
#
FUNCTION(TRIBITS_DETERMINE_IF_CURRENT_PACKAGE_NEEDS_REBUILT)

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(TIMER_START_SECONDS)
  ENDIF()

   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "CURRENT_PACKAGE_OUT_OF_DATE_OUT"
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     ${ARGN}
     )

  # Get pass through print level options 
  SET(SHOW_MOST_RECENT_FILES_ARGS)
  IF (PARSE_SHOW_MOST_RECENT_FILES)
    APPEND_SET(SHOW_MOST_RECENT_FILES_ARGS SHOW_MOST_RECENT_FILES)
  ENDIF()
  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    APPEND_SET(SHOW_MOST_RECENT_FILES_ARGS SHOW_OVERALL_MOST_RECENT_FILES)
  ENDIF()
  #PRINT_VAR(SHOW_MOST_RECENT_FILES_ARGS)

  IF (PARSE_SHOW_MOST_RECENT_FILES)
    SET(PARSE_SHOW_OVERALL_MOST_RECENT_FILES TRUE)
  ENDIF()

  #
  # B) Get the list of enabled upstream packages
  #

  # Only search parent packages to cut down on dirs searched
  SET(ENABLED_UPSTREAM_PACKAGES)
  SET(CURRENT_PARENT_PACKAGE)
  FOREACH(UPSTREAM_SE_PACKAGE ${${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES})
    # Assume we will append
    SET(APPEND_PACKAGE ${UPSTREAM_SE_PACKAGE})
    # If is a subpackage we only append the parent packages
    SET(PARENT_PACKAGE ${${UPSTREAM_SE_PACKAGE}_PARENT_PACKAGE})
    IF (PARENT_PACKAGE)
      SET(APPEND_PACKAGE ${PARENT_PACKAGE})
    ENDIF()
    # Append
    APPEND_SET(ENABLED_UPSTREAM_PACKAGES ${APPEND_PACKAGE})
  ENDFOREACH()
  LIST(REMOVE_DUPLICATES ENABLED_UPSTREAM_PACKAGES)
  #PRINT_VAR(ENABLED_UPSTREAM_PACKAGES)

  #
  # C) Determine the most recent files on the upstream SE packages
  #

  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    MESSAGE("\nDetermining most recent source file in upstream SE packages"
      " from ${PACKAGE_NAME}:")  
  ENDIF()
  SET(UPSTREAM_SOURCE_BASE_DIRS)
  FOREACH(UPSTREAM_PACKAGE ${ENABLED_UPSTREAM_PACKAGES})
    APPEND_SET(UPSTREAM_SOURCE_BASE_DIRS "${${UPSTREAM_PACKAGE}_SOURCE_DIR}") 
  ENDFOREACH()
  TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP(
    SOURCE_BASE_DIRS ${UPSTREAM_SOURCE_BASE_DIRS}
    SOURCE_BASE_BASE_DIR "${PROJECT_SOURCE_DIR}"
    ${SHOW_MOST_RECENT_FILES_ARGS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT MOST_RECENT_UPSTREAM_SOURCE_FILEPATH
    )
  #PRINT_VAR(MOST_RECENT_UPSTREAM_SOURCE_FILEPATH)

  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    MESSAGE("\nDetermining most recent binary file in upstream SE packages"
      " from ${PACKAGE_NAME}:")  
  ENDIF()
  SET(UPSTREAM_BINARY_BASE_DIRS)
  FOREACH(UPSTREAM_PACKAGE ${ENABLED_UPSTREAM_PACKAGES})
    APPEND_SET(UPSTREAM_BINARY_BASE_DIRS "${${UPSTREAM_PACKAGE}_BINARY_DIR}") 
  ENDFOREACH()
  TRIBITS_FIND_MOST_RECENT_BINARY_FILE_TIMESTAMP(
    BINARY_BASE_DIRS ${UPSTREAM_BINARY_BASE_DIRS}
    BINARY_BASE_BASE_DIR "${PROJECT_BINARY_DIR}"
    ${SHOW_MOST_RECENT_FILES_ARGS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_UPSTREAM_BINARY_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT MOST_RECENT_UPSTREAM_BINARY_FILEPATH
    )
  #PRINT_VAR(MOST_RECENT_UPSTREAM_BINARY_FILEPATH)

  #
  # D) Determine the most recent files for the current package
  #

  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    MESSAGE("\nDetermining most recent source file for current"
      " package ${PACKAGE_NAME}:")  
  ENDIF()
  TRIBITS_FIND_MOST_RECENT_SOURCE_FILE_TIMESTAMP(
    SOURCE_BASE_DIRS ${${PACKAGE_NAME}_SOURCE_DIR}
    SOURCE_BASE_BASE_DIR "${PROJECT_SOURCE_DIR}"
    ${SHOW_MOST_RECENT_FILES_ARGS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_THIS_PACKAGE_SOURCE_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT MOST_RECENT_THIS_SOURCE_FILEPATH
    )

  IF (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    MESSAGE("\nDetermining most recent binary file for current"
      " package ${PACKAGE_NAME}:")  
  ENDIF()
  TRIBITS_FIND_MOST_RECENT_BINARY_FILE_TIMESTAMP(
    BINARY_BASE_DIRS  ${${PACKAGE_NAME}_BINARY_DIR}
    BINARY_BASE_BASE_DIR "${PROJECT_BINARY_DIR}"
    ${SHOW_MOST_RECENT_FILES_ARGS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_THIS_PACKAGE_BINARY_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT  MOST_RECENT_THIS_BINARY_FILEPATH
    )

  #
  # E) Compare most recent file time stamps to determine if a rebuild is needed
  # 

  SET(CURRENT_PACKAGE_OUT_OF_DATE_OUT FALSE)

  MESSAGE("\nComparing timestamps of recently updated files:")

  IF (MOST_RECENT_THIS_BINARY_FILEPATH)

    TRIBITS_UPDATE_PACKAGE_OUT_OF_DATE(
      "upstream SE package source" ${MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP}
         "${MOST_RECENT_UPSTREAM_SOURCE_FILEPATH}"
      ${MOST_RECENT_THIS_PACKAGE_BINARY_TIMESTAMP} "${MOST_RECENT_THIS_BINARY_FILEPATH}"
      CURRENT_PACKAGE_OUT_OF_DATE_OUT )
  
    TRIBITS_UPDATE_PACKAGE_OUT_OF_DATE(
      "upstream SE package binary" ${MOST_RECENT_UPSTREAM_BINARY_TIMESTAMP}
         "${MOST_RECENT_UPSTREAM_BINARY_FILEPATH}"
      ${MOST_RECENT_THIS_PACKAGE_BINARY_TIMESTAMP} "${MOST_RECENT_THIS_BINARY_FILEPATH}"
      CURRENT_PACKAGE_OUT_OF_DATE_OUT )
  
    TRIBITS_UPDATE_PACKAGE_OUT_OF_DATE(
      "this package's source" ${MOST_RECENT_THIS_PACKAGE_SOURCE_TIMESTAMP}
         "${MOST_RECENT_THIS_SOURCE_FILEPATH}"
      ${MOST_RECENT_THIS_PACKAGE_BINARY_TIMESTAMP} "${MOST_RECENT_THIS_BINARY_FILEPATH}"
      CURRENT_PACKAGE_OUT_OF_DATE_OUT )
  
    IF (NOT CURRENT_PACKAGE_OUT_OF_DATE_OUT)
      MESSAGE("-- This package's most recent binary file"
        " ${MOST_RECENT_THIS_BINARY_FILEPATH}"
        " is more recent than its upstream SE package source or binary files"
        " or this package's source files!")
    ENDIF()

  ELSE()

    MESSAGE("-- This package has no unfiltered binary files so consider out of date!")

  ENDIF()

  SET(${PARSE_CURRENT_PACKAGE_OUT_OF_DATE_OUT} ${CURRENT_PACKAGE_OUT_OF_DATE_OUT}
    PARENT_SCOPE)

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(TIMER_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${TIMER_START_SECONDS} ${TIMER_STOP_SECONDS}
      "\nTotal time to check for most recent modified files")
  ENDIF()
   
ENDFUNCTION()


#
# Utility functions
#


FUNCTION(TRIBITS_UPDATE_PACKAGE_OUT_OF_DATE
  DEPENDENCY_TYPE_STRING  DEP_FILE_TIMESTAMP  DEP_FILEPATH
  THIS_BINARY_FILE_TIMESTAMP  THIS_BINARY_FILEPATH
  CURRENT_PACKAGE_IS_OUT_OF_DATE_INOUT
  )
  IF ("${DEP_FILE_TIMESTAMP}" GREATER "${THIS_BINARY_FILE_TIMESTAMP}")
    MESSAGE("-- The ${DEPENDENCY_TYPE_STRING} file ${DEP_FILEPATH} is more recent than"
      " this package's binary file ${THIS_BINARY_FILEPATH}!")
    SET(${CURRENT_PACKAGE_IS_OUT_OF_DATE_INOUT} TRUE PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()


FUNCTION(TRIBITS_GET_HUMAN_READABLE_FILE_AND_TIMESTAMP
  BASE_DIR   CURRENT_FILEPATH
  HUMAN_READABLE_FILE_AND_TIMESTAMP_OUT
  )
  EXECUTE_PROCESS(
    WORKING_DIRECTORY "${BASE_DIR}"
    COMMAND ls --full-time "${CURRENT_FILEPATH}"
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE  HUMAN_READABLE_FILE_AND_TIMESTAMP
    )
  SET(${HUMAN_READABLE_FILE_AND_TIMESTAMP_OUT}
    ${HUMAN_READABLE_FILE_AND_TIMESTAMP}
    PARENT_SCOPE)
ENDFUNCTION()





# LocalWords:  ENDFOREACH subpackage subpackages TriBITS timestamp
