# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(TribitsConfigureTiming)

include(CMakeParseArguments)


# @FUNCTION: tribits_find_most_recent_file_timestamp()
#
# Find the most modified file in a set of base directories and return its
# timestamp.
#
# Usage::
#
#   tribits_find_most_recent_file_timestamp(
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
#     $ find . -type f -printf '%T@ %p\n' \
#         | grep -v "<re0>" | grep -v "<re1>" | ... \
#         | sort -n | tail -1
#
# to return the most recent file in each listed directory ``<dir0>``,
# ``<dir1>``, etc.  It then determines the most recently modified file over
# all of the directories and prints and returns in the variables
# ``<mostRecentTimestamp>``, ``<mostRecentFilepathBaseDir>``, and
# ``<mostRecentRelativeFilePath>``.
#
function(tribits_find_most_recent_file_timestamp)

    if (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
      message("\nSearching for most modified files in base dirs:")
    endif()

  #
  # A) Parse the input arguments
  #

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     #one_value_keywords
     ""
     #mulit_value_keywords
     "BASE_DIRS;BASE_BASE_DIR;EXCLUDE_REGEXES;MOST_RECENT_TIMESTAMP_OUT;MOST_RECENT_FILEPATH_BASE_DIR_OUT;MOST_RECENT_RELATIVE_FILEPATH_OUT"
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  if (PARSE_SHOW_MOST_RECENT_FILES)
    set(PARSE_SHOW_OVERALL_MOST_RECENT_FILE ON)
  endif()

  #
  # B) Loop over each directory and find the most modified file
  #

  set(OVERALL_MOST_RECENT_TIMESTAMP "0000000000.0000000000")
  set(OVERALL_MOST_RECENT_FILEPATH "")
  set(OVERALL_MOST_RECENT_FILEPATH_DIR "")
  set(OVERALL_MOST_RECENT_RELATEIVE_FILEPATH_DIR "")
  set(OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE "")

  foreach(BASE_DIR ${PARSE_BASE_DIRS})

    if (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
      message("\nSearching '${BASE_DIR}' ...")
    endif()

    if (IS_DIRECTORY "${BASE_DIR}")

      # Build up commands for grep -v
      set(GREP_V_COMMANDS)
      foreach(EXCLUDE_REGEX ${PARSE_EXCLUDE_REGEXES})
        append_set(GREP_V_COMMANDS COMMAND grep -v "${EXCLUDE_REGEX}")
      endforeach()

      # Get the time stamp and the file name of the most recently modified file
      # in current directory.
      execute_process(
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

      if (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
        print_var(MOST_RECENT_TIMESTAMP_AND_FILE)
      endif()

      if (MOST_RECENT_TIMESTAMP_AND_FILE)

        split("${MOST_RECENT_TIMESTAMP_AND_FILE}" " "
          MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT)

        # Get the time stamp part
        list(GET MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT 0
          CURRENT_TIMESTAMP)

        # Get the relative file path
        list(GET MOST_RECENT_TIMESTAMP_AND_FILE_SPLIT 1
          CURRENT_FILEPATH)

        # Get the directory relative to the base base dir
        string(REPLACE "${PARSE_BASE_BASE_DIR}/" "./"  RELATIVE_FILEPATH_DIR
          "${BASE_DIR}")

        if (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
          print_var(CURRENT_TIMESTAMP)
          print_var(CURRENT_FILEPATH)
        endif()

        if (PARSE_SHOW_MOST_RECENT_FILES)
          tribits_get_human_readable_file_and_timestamp(
            "${BASE_DIR}"  "${CURRENT_FILEPATH}"
            HUMAN_READABLE_FILE_AND_TIMESTAMP
            )
          message("-- " "Most recent file in ${RELATIVE_FILEPATH_DIR}/"
            " is ${CURRENT_FILEPATH}\n"
            "    ${HUMAN_READABLE_FILE_AND_TIMESTAMP}")
        endif()

        if ("${CURRENT_TIMESTAMP}" GREATER "${OVERALL_MOST_RECENT_TIMESTAMP}")
          if (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
            message("    New most recent file path!")
          endif()
          set(OVERALL_MOST_RECENT_TIMESTAMP "${CURRENT_TIMESTAMP}")
          set(OVERALL_MOST_RECENT_RELATIVE_FILEPATH "${CURRENT_FILEPATH}")
          set(OVERALL_MOST_RECENT_FILEPATH_DIR "${BASE_DIR}")
          set(OVERALL_MOST_RECENT_RELATEIVE_FILEPATH_DIR "${RELATIVE_FILEPATH_DIR}")
          if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
            set(OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE
              "${HUMAN_READABLE_FILE_AND_TIMESTAMP}")
          endif()
        endif()

      endif()

    else()

      if (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
        message("Directory does not exist, skipping ...")
      endif()

    endif()

  endforeach()

  if (TRIBITS_FIND_MOST_RECENT_FILE_TIMESTAMP_DEBUG_DUMP)
    print_var(OVERALL_MOST_RECENT_TIMESTAMP)
    print_var(OVERALL_MOST_RECENT_RELATIVE_FILEPATH)
  endif()

  if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    if (OVERALL_MOST_RECENT_FILEPATH_DIR)
      if (NOT OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE)
        tribits_get_human_readable_file_and_timestamp(
          "${OVERALL_MOST_RECENT_FILEPATH_DIR}"
          "${OVERALL_MOST_RECENT_RELATIVE_FILEPATH}"
          OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE
        )
      endif()
      message("-- " "Overall most recent modified file is in"
        " ${OVERALL_MOST_RECENT_RELATEIVE_FILEPATH_DIR}/ and is ${OVERALL_MOST_RECENT_RELATIVE_FILEPATH}\n"
        "    ${OVERALL_MOST_RECENT_FILEPATH_TIMESTAMP_HUMAN_READABLE}")
    else()
      message("-- There are no unfiltered files!")
    endif()
  endif()

  set(${PARSE_MOST_RECENT_TIMESTAMP_OUT}
    ${OVERALL_MOST_RECENT_TIMESTAMP}
    PARENT_SCOPE)

  if (PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT)
    set(${PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT}
      ${OVERALL_MOST_RECENT_FILEPATH_DIR}
      PARENT_SCOPE)
  endif()

  if (PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT)
    set(${PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT}
      ${OVERALL_MOST_RECENT_RELATIVE_FILEPATH}
      PARENT_SCOPE )
  endif()

endfunction()


# @FUNCTION: tribits_find_most_recent_source_file_timestamp()
#
# Find the most modified source file in a set of base directories and return
# its timestamp.
#
# Usage::
#
#   tribits_find_most_recent_source_file_timestamp(
#     SOURCE_BASE_DIRS <dir0> <dir1> ...
#     [SOURCE_BASE_BASE_DIR <dir>]
#     [SHOW_MOST_RECENT_FILES]
#     [SHOW_OVERALL_MOST_RECENT_FILE]
#     [MOST_RECENT_TIMESTAMP_OUT  <mostRecentTimestamp>]
#     [MOST_RECENT_FILEPATH_BASE_DIR_OUT <mostRecentFilepathBaseDir>]
#     [MOST_RECENT_RELATIVE_FILEPATH_OUT <mostRecentRelativeFilePath>]
#     )
#
# This function just calls `tribits_find_most_recent_file_timestamp()`_
# passing in a set of basic exclude regexes like ``[.]git/``, ``[.]svn/``,
# etc.  These types of version control files can not possibly directly impact
# the source code.
#
function(tribits_find_most_recent_source_file_timestamp)

  #
  # A) Parse the input arguments
  #

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     #one_value_keywords
     ""
     #multi_value_keywords
     "SOURCE_BASE_DIRS;SOURCE_BASE_BASE_DIR;MOST_RECENT_TIMESTAMP_OUT;MOST_RECENT_FILEPATH_BASE_DIR_OUT;MOST_RECENT_RELATIVE_FILEPATH_OUT"
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  #
  # B) Call the function tribits_find_most_recent_file_timestamp()
  #

  set(FILTER_OUT_SOURCE_FILE_REGEXS
    "/[.]git/"
    )

  set(VARIABLE_ARGS)
  if (PARSE_SHOW_MOST_RECENT_FILES)
    append_set(VARIABLE_ARGS SHOW_MOST_RECENT_FILES)
  endif()
  if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    append_set(VARIABLE_ARGS SHOW_OVERALL_MOST_RECENT_FILES)
  endif()

  #print_var(VARIABLE_ARGS)
  tribits_find_most_recent_file_timestamp(
    BASE_DIRS ${PARSE_SOURCE_BASE_DIRS}
    BASE_BASE_DIR ${PARSE_SOURCE_BASE_BASE_DIR}
    EXCLUDE_REGEXES ${FILTER_OUT_SOURCE_FILE_REGEXS}
    MOST_RECENT_TIMESTAMP_OUT MOST_RECENT_TIMESTAMP
    MOST_RECENT_FILEPATH_BASE_DIR_OUT  MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT  MOST_RECENT_RELATIVE_FILEPATH
    ${VARIABLE_ARGS}
    )
  #print_var(MOST_RECENT_TIMESTAMP)

  set(${PARSE_MOST_RECENT_TIMESTAMP_OUT} ${MOST_RECENT_TIMESTAMP}
    PARENT_SCOPE)

  if (PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT)
    set(${PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT} ${MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP}
      PARENT_SCOPE)
  endif()

  if (PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT)
    set(${PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT} ${MOST_RECENT_RELATIVE_FILEPATH}
      PARENT_SCOPE )
  endif()

endfunction()


# @FUNCTION: tribits_find_most_recent_binary_file_timestamp()
#
# Find the most modified binary file in a set of base directories and return
# its timestamp.
#
# Usage::
#
#   tribits_find_most_recent_binary_file_timestamp(
#     BINARY_BASE_DIRS <dir0> <dir1> ...
#     [BINARY_BASE_BASE_DIR <dir>]
#     [MOST_RECENT_TIMESTAMP_OUT  <mostRecentTimestamp>]
#     [MOST_RECENT_FILEPATH_BASE_DIR_OUT <mostRecentFilepathBaseDir>]
#     [MOST_RECENT_RELATIVE_FILEPATH_OUT <mostRecentRelativeFilePath>]
#     [SHOW_MOST_RECENT_FILES]
#     [SHOW_OVERALL_MOST_RECENT_FILE]
#     )
#
# This function just calls `tribits_find_most_recent_file_timestamp()`_
# passing in a set of basic exclude regexes like ``CMakeFiles/``,
# ``[.]cmake$``, and ``/Makefile$``, etc.  These types of files usually don't
# impact the build of downstream software in CMake projects.
#
function(tribits_find_most_recent_binary_file_timestamp)

  #
  # A) Parse the input arguments
  #

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     #one_value_keywords
     ""
     #multi_value_keywords
     "BINARY_BASE_DIRS;BINARY_BASE_BASE_DIR;MOST_RECENT_TIMESTAMP_OUT;MOST_RECENT_FILEPATH_BASE_DIR_OUT;MOST_RECENT_RELATIVE_FILEPATH_OUT"
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  #
  # B) Define filters for binary files we know are not significant
  #

  set(FILTER_OUT_BINARY_FILE_REGEXS
    "CMakeFiles/" "[.]cmake$" "/Makefile$"
    )

  #
  # C) Call the function tribits_find_most_recent_file_timestamp()
  #

  set(VARIABLE_ARGS)
  if (PARSE_SHOW_MOST_RECENT_FILES)
    append_set(VARIABLE_ARGS SHOW_MOST_RECENT_FILES)
  endif()
  if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    append_set(VARIABLE_ARGS SHOW_OVERALL_MOST_RECENT_FILES)
  endif()

  #print_var(VARIABLE_ARGS)
  tribits_find_most_recent_file_timestamp(
    BASE_DIRS ${PARSE_BINARY_BASE_DIRS}
    BASE_BASE_DIR ${PARSE_BINARY_BASE_BASE_DIR}
    EXCLUDE_REGEXES ${FILTER_OUT_BINARY_FILE_REGEXS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT  MOST_RECENT_RELATIVE_FILEPATH
    ${VARIABLE_ARGS}
    )
  #print_var(MOST_RECENT_TIMESTAMP)

  set(${PARSE_MOST_RECENT_TIMESTAMP_OUT} ${MOST_RECENT_TIMESTAMP}
    PARENT_SCOPE)

  if (PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT)
    set(${PARSE_MOST_RECENT_FILEPATH_BASE_DIR_OUT} ${MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP}
      PARENT_SCOPE)
  endif()

  if (PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT)
    set(${PARSE_MOST_RECENT_RELATIVE_FILEPATH_OUT} ${MOST_RECENT_RELATIVE_FILEPATH}
      PARENT_SCOPE )
  endif()

endfunction()


# @FUNCTION: tribits_determine_if_current_package_needs_rebuilt()
#
# Determine at configure time if any of the upstream dependencies for a
# package require the current package to be rebuilt.
#
# Usage::
#
#   tribits_determine_if_current_package_needs_rebuilt(
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
#     individual directories for each category (i.e. one for upstream package
#     source dirs, one for upstream package binary dirs, one for the package's
#     source dir, and one for the package's own binary dir) is printed to
#     STDOUT.
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
# * `tribits_find_most_recent_source_file_timestamp()`_
# * `tribits_find_most_recent_binary_file_timestamp()`_
#
# to determine the most recent modified files in the upstream TriBITS
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
function(tribits_determine_if_current_package_needs_rebuilt)

  tribits_config_code_start_timer(TIMER_START_SECONDS)

  #
  # A) Parse the input arguments
  #

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     "SHOW_MOST_RECENT_FILES;SHOW_OVERALL_MOST_RECENT_FILES"
     #one_value_keywords
     ""
     #mulit_value_keywords
     "CURRENT_PACKAGE_OUT_OF_DATE_OUT"
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  # Get pass through print level options
  set(SHOW_MOST_RECENT_FILES_ARGS)
  if (PARSE_SHOW_MOST_RECENT_FILES)
    append_set(SHOW_MOST_RECENT_FILES_ARGS SHOW_MOST_RECENT_FILES)
  endif()
  if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    append_set(SHOW_MOST_RECENT_FILES_ARGS SHOW_OVERALL_MOST_RECENT_FILES)
  endif()
  #print_var(SHOW_MOST_RECENT_FILES_ARGS)

  if (PARSE_SHOW_MOST_RECENT_FILES)
    set(PARSE_SHOW_OVERALL_MOST_RECENT_FILES TRUE)
  endif()

  #
  # B) Get the list of enabled upstream packages
  #

  # Only search parent packages to cut down on dirs searched
  set(ENABLED_UPSTREAM_PACKAGES)
  set(CURRENT_PARENT_PACKAGE)
  foreach(upstreamPackage ${${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES})
    # Assume we will append
    set(APPEND_PACKAGE ${upstreamPackage})
    # If is a subpackage we only append the parent packages
    set(PARENT_PACKAGE ${${upstreamPackage}_PARENT_PACKAGE})
    if (PARENT_PACKAGE)
      set(APPEND_PACKAGE ${PARENT_PACKAGE})
    endif()
    # Append
    append_set(ENABLED_UPSTREAM_PACKAGES ${APPEND_PACKAGE})
  endforeach()
  list(REMOVE_DUPLICATES ENABLED_UPSTREAM_PACKAGES)
  #print_var(ENABLED_UPSTREAM_PACKAGES)

  #
  # C) Determine the most recent files on the upstream packages
  #

  if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    message("\nDetermining most recent source file in upstream packages"
      " from ${PACKAGE_NAME}:")
  endif()
  set(UPSTREAM_SOURCE_BASE_DIRS)
  foreach(UPSTREAM_PACKAGE ${ENABLED_UPSTREAM_PACKAGES})
    append_set(UPSTREAM_SOURCE_BASE_DIRS "${${UPSTREAM_PACKAGE}_SOURCE_DIR}")
  endforeach()
  tribits_find_most_recent_source_file_timestamp(
    SOURCE_BASE_DIRS ${UPSTREAM_SOURCE_BASE_DIRS}
    SOURCE_BASE_BASE_DIR "${PROJECT_SOURCE_DIR}"
    ${SHOW_MOST_RECENT_FILES_ARGS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT MOST_RECENT_UPSTREAM_SOURCE_FILEPATH
    )
  #print_var(MOST_RECENT_UPSTREAM_SOURCE_FILEPATH)

  if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    message("\nDetermining most recent binary file in upstream packages"
      " from ${PACKAGE_NAME}:")
  endif()
  set(UPSTREAM_BINARY_BASE_DIRS)
  foreach(UPSTREAM_PACKAGE ${ENABLED_UPSTREAM_PACKAGES})
    append_set(UPSTREAM_BINARY_BASE_DIRS "${${UPSTREAM_PACKAGE}_BINARY_DIR}")
  endforeach()
  tribits_find_most_recent_binary_file_timestamp(
    BINARY_BASE_DIRS ${UPSTREAM_BINARY_BASE_DIRS}
    BINARY_BASE_BASE_DIR "${PROJECT_BINARY_DIR}"
    ${SHOW_MOST_RECENT_FILES_ARGS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_UPSTREAM_BINARY_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT MOST_RECENT_UPSTREAM_BINARY_FILEPATH
    )
  #print_var(MOST_RECENT_UPSTREAM_BINARY_FILEPATH)

  #
  # D) Determine the most recent files for the current package
  #

  if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    message("\nDetermining most recent source file for current"
      " package ${PACKAGE_NAME}:")
  endif()
  tribits_find_most_recent_source_file_timestamp(
    SOURCE_BASE_DIRS ${${PACKAGE_NAME}_SOURCE_DIR}
    SOURCE_BASE_BASE_DIR "${PROJECT_SOURCE_DIR}"
    ${SHOW_MOST_RECENT_FILES_ARGS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_THIS_PACKAGE_SOURCE_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT MOST_RECENT_THIS_SOURCE_FILEPATH
    )

  if (PARSE_SHOW_OVERALL_MOST_RECENT_FILES)
    message("\nDetermining most recent binary file for current"
      " package ${PACKAGE_NAME}:")
  endif()
  tribits_find_most_recent_binary_file_timestamp(
    BINARY_BASE_DIRS  ${${PACKAGE_NAME}_BINARY_DIR}
    BINARY_BASE_BASE_DIR "${PROJECT_BINARY_DIR}"
    ${SHOW_MOST_RECENT_FILES_ARGS}
    MOST_RECENT_TIMESTAMP_OUT  MOST_RECENT_THIS_PACKAGE_BINARY_TIMESTAMP
    MOST_RECENT_RELATIVE_FILEPATH_OUT  MOST_RECENT_THIS_BINARY_FILEPATH
    )

  #
  # E) Compare most recent file time stamps to determine if a rebuild is needed
  #

  set(CURRENT_PACKAGE_OUT_OF_DATE_OUT FALSE)

  message("\nComparing timestamps of recently updated files:")

  if (MOST_RECENT_THIS_BINARY_FILEPATH)

    tribits_update_package_out_of_date(
      "upstream package source" ${MOST_RECENT_UPSTREAM_SOURCE_TIMESTAMP}
         "${MOST_RECENT_UPSTREAM_SOURCE_FILEPATH}"
      ${MOST_RECENT_THIS_PACKAGE_BINARY_TIMESTAMP} "${MOST_RECENT_THIS_BINARY_FILEPATH}"
      CURRENT_PACKAGE_OUT_OF_DATE_OUT )

    tribits_update_package_out_of_date(
      "upstream package binary" ${MOST_RECENT_UPSTREAM_BINARY_TIMESTAMP}
         "${MOST_RECENT_UPSTREAM_BINARY_FILEPATH}"
      ${MOST_RECENT_THIS_PACKAGE_BINARY_TIMESTAMP} "${MOST_RECENT_THIS_BINARY_FILEPATH}"
      CURRENT_PACKAGE_OUT_OF_DATE_OUT )

    tribits_update_package_out_of_date(
      "this package's source" ${MOST_RECENT_THIS_PACKAGE_SOURCE_TIMESTAMP}
         "${MOST_RECENT_THIS_SOURCE_FILEPATH}"
      ${MOST_RECENT_THIS_PACKAGE_BINARY_TIMESTAMP} "${MOST_RECENT_THIS_BINARY_FILEPATH}"
      CURRENT_PACKAGE_OUT_OF_DATE_OUT )

    if (NOT CURRENT_PACKAGE_OUT_OF_DATE_OUT)
      message("-- This package's most recent binary file"
        " ${MOST_RECENT_THIS_BINARY_FILEPATH}"
        " is more recent than its upstream package source or binary files"
        " or this package's source files!")
    endif()

  else()

    message("-- This package has no unfiltered binary files so consider out of date!")

  endif()

  set(${PARSE_CURRENT_PACKAGE_OUT_OF_DATE_OUT} ${CURRENT_PACKAGE_OUT_OF_DATE_OUT}
    PARENT_SCOPE)

  tribits_config_code_stop_timer(TIMER_START_SECONDS
    "\nTotal time to check for most recent modified files")

endfunction()


#
# Utility functions
#


function(tribits_update_package_out_of_date
  DEPENDENCY_TYPE_STRING  DEP_FILE_TIMESTAMP  DEP_FILEPATH
  THIS_BINARY_FILE_TIMESTAMP  THIS_BINARY_FILEPATH
  CURRENT_PACKAGE_IS_OUT_OF_DATE_INOUT
  )
  if ("${DEP_FILE_TIMESTAMP}" GREATER "${THIS_BINARY_FILE_TIMESTAMP}")
    message("-- The ${DEPENDENCY_TYPE_STRING} file ${DEP_FILEPATH} is more recent than"
      " this package's binary file ${THIS_BINARY_FILEPATH}!")
    set(${CURRENT_PACKAGE_IS_OUT_OF_DATE_INOUT} TRUE PARENT_SCOPE)
  endif()
endfunction()


function(tribits_get_human_readable_file_and_timestamp
  BASE_DIR   CURRENT_FILEPATH
  HUMAN_READABLE_FILE_AND_TIMESTAMP_OUT
  )
  execute_process(
    WORKING_DIRECTORY "${BASE_DIR}"
    COMMAND ls --full-time "${CURRENT_FILEPATH}"
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE  HUMAN_READABLE_FILE_AND_TIMESTAMP
    )
  set(${HUMAN_READABLE_FILE_AND_TIMESTAMP_OUT}
    ${HUMAN_READABLE_FILE_AND_TIMESTAMP}
    PARENT_SCOPE)
endfunction()





# LocalWords:  ENDFOREACH subpackage subpackages TriBITS timestamp
