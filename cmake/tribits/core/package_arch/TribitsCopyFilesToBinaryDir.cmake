# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(CMakeParseArguments)


# @FUNCTION: tribits_copy_files_to_binary_dir()
#
# Function that copies a list of files from a source directory to a
# destination directory at configure time, typically so that it can be used in
# one or more tests.
#
# Usage::
#
#   tribits_copy_files_to_binary_dir(
#     <targetName>
#     [SOURCE_FILES <file1> <file2> ...]
#     [SOURCE_DIR <sourceDir>]
#     [DEST_FILES <dfile1> <dfile2> ...]
#     [DEST_DIR <destDir>]
#     [TARGETDEPS <targDep1> <targDep2> ...]
#     [EXEDEPS <exeDep1> <exeDep2> ...]
#     [NOEXEPREFIX]
#     [CATEGORIES <category1>  <category2> ...]
#     )
#
# This sets up all of the custom CMake commands and targets to ensure that the
# files in the destination directory are always up to date just by building
# the ``ALL`` target.
#
# **NOTE:** The target name ``<targetName>`` must be unique from all other
# targets in the same TriBITS Package.  Otherwise, one will get a configure
# failure complaining that a target name has already been defined.  Therefore,
# be sure to pick long and unique target names!
#
# This function has a few valid calling modes:
#
# **1) Source files and destination files have the same name**::
#
#   tribits_copy_files_to_binary_dir(
#     <targetName>
#     SOURCE_FILES <file1> <file2> ...
#     [SOURCE_DIR <sourceDir>]
#     [DEST_DIR <destDir>]
#     [TARGETDEPS <targDep1> <targDep2> ...]
#     [EXEDEPS <exeDep1> <exeDep2> ...]
#     [NOEXEPREFIX]
#     [CATEGORIES <category1>  <category2> ...]
#     )
#
# In this case, the names of the source files and the destination files
# are the same but just live in different directories.
#
# **2) Source files have a prefix different from the destination files**::
#
#   tribits_copy_files_to_binary_dir(
#     <targetName>
#     DEST_FILES <file1> <file2> ...
#     SOURCE_PREFIX <srcPrefix>
#     [SOURCE_DIR <sourceDir>]
#     [DEST_DIR <destDir>]
#     [EXEDEPS <exeDep1> <exeDep2> ...]
#     [NOEXEPREFIX]
#     [CATEGORIES <category1>  <category2> ...]
#     )
#
# In this case, the source files have the same basic name as the destination
# files except they have the prefix ``<srcPrefix>`` prepended to the name.
#
# **3) Source files and destination files have completely different names**::
#
#   tribits_copy_files_to_binary_dir(
#     <targetName>
#     SOURCE_FILES <sfile1> <sfile2> ...
#     [SOURCE_DIR <sourceDir>]
#     DEST_FILES <dfile1> <dfile2> ...
#     [DEST_DIR <destDir>]
#     [EXEDEPS <exeDep1> <exeDep2> ...]
#     [NOEXEPREFIX]
#     [CATEGORIES <category1>  <category2> ...]
#     )
#
# In this case, the source files and destination files have completely
# different prefixes.
#
# The individual arguments are:
#
#   ``SOURCE_FILES <file1> <file2> ...``
#
#     Listing of the source files relative to the source directory given by
#     the argument ``SOURCE_DIR <sourceDir>``.  If omitted, this list will be
#     the same as ``DEST_FILES`` with the argument ``SOURCE_PREFIX
#     <srcPrefix>`` appended.
#
#   ``SOURCE_DIR <sourceDir>``
#
#     Optional argument that gives the (absolute) base directory for all of
#     the source files.  If omitted, this takes the default value of
#     ``${CMAKE_CURRENT_SOURCE_DIR}``.
#
#   ``DEST_FILES <file1> <file2> ...``
#
#     Listing of the destination files relative to the destination directory
#     given by the argument ``DEST_DIR <destDir>``. If omitted, this list will
#     be the same as given by the ``SOURCE_FILES`` list.
#
#   ``DEST_DIR <destDir>``
#
#     Optional argument that gives the (absolute) base directory for all of
#     the destination files.  If omitted, this takes the default value of
#     ``${CMAKE_CURRENT_BINARY_DIR}``
#
#   ``TARGETDEPS <targDep1> <targDep2> ...``
#
#     Listing of general CMake targets that these files will be added as
#     dependencies to.  This results in the copies to be performed when any of
#     the targets ``<targDepi>`` are built.
#
#   ``EXEDEPS <exeDep1> <exeDep2> ...``
#
#     Listing of executable targets that these files will be added as
#     dependencies to.  By default, the prefix ``${PACKAGE_NAME}_`` will is
#     appended to the names of the targets.  This ensures that if the
#     executable target is built that these files will also be copied as well.
#
#   ``NOEXEPREFIX``
#
#     Option that determines if the prefix ``${PACKAGE_NAME}_`` will be
#     appended to the arguments in the ``EXEDEPS`` list.
#
function(tribits_copy_files_to_binary_dir TARGET_NAME)

  #
  # A) Parse input arguments
  #
  cmake_parse_arguments(
    #prefix
    PARSE
    #options
    "NOEXEPREFIX"
    #one_value_keywords
    ""
    #multi_value_keywords
    "SOURCE_DIR;SOURCE_FILES;SOURCE_PREFIX;DEST_DIR;DEST_FILES;EXEDEPS;TARGETDEPS;CATEGORIES"
    ${ARGN}
  )

  tribits_check_for_unparsed_arguments()

  set(ADD_THE_TEST FALSE)
  tribits_add_test_process_categories(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    return()
  endif()

  if (NOT PARSE_SOURCE_DIR)
    set(PARSE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  if (NOT PARSE_DEST_DIR)
    set(PARSE_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  set(RENAME 1)

  if (PARSE_SOURCE_FILES AND NOT PARSE_DEST_FILES)
    set(PARSE_DEST_FILES ${PARSE_SOURCE_FILES})
    set(RENAME 0)
  endif()

  if (PARSE_DEST_FILES AND NOT PARSE_SOURCE_FILES)
    set(PARSE_SOURCE_FILES ${PARSE_DEST_FILES})
    set(RENAME 0)
  endif()

  #
  # B) Validate arguments
  #

  list(LENGTH PARSE_SOURCE_DIR PARSE_SOURCE_DIR_LEN)
  list(LENGTH PARSE_SOURCE_FILES PARSE_SOURCE_FILES_LEN)
  list(LENGTH PARSE_SOURCE_PREFIX PARSE_SOURCE_PREFIX_LEN)
  list(LENGTH PARSE_DEST_DIR PARSE_DEST_DIR_LEN)
  list(LENGTH PARSE_DEST_FILES PARSE_DEST_FILES_LEN)

  if (PARSE_SOURCE_DIR_LEN GREATER 1)
    message(SEND_ERROR "Error, there can only be 0 or one SOURCE_DIR arguments!")
  endif()

  if (PARSE_DEST_DIR_LEN GREATER 1)
    message(SEND_ERROR "Error, there can only be 0 or one DEST_DIR arguments!")
  endif()

  if (PARSE_SOURCE_PREFIX_LEN GREATER 1)
    message(SEND_ERROR "Error, If SOURCE_PREFIX can only take one argument!")
  endif()

  if (PARSE_SOURCE_FILES_LEN EQUAL 0)
    message(SEND_ERROR "Error, there are no source files listed!")
  endif()

  if (PARSE_DEST_FILES_LEN EQUAL 0)
    message(SEND_ERROR "Error, there are no destination files listed!")
  endif()

  if (NOT PARSE_SOURCE_FILES_LEN EQUAL ${PARSE_DEST_FILES_LEN})
    message(SEND_ERROR "Error, there are not the same number of source files ${PARSE_SOURCE_FILES_LEN} and dest files ${PARSE_DEST_FILES_LEN}!")
  endif()

  #
  # C) Build the list of command and dependencies
  #

  set(DEST_FILES_LIST)
  if(RENAME)
    math(EXPR FILES_IDX_END "${PARSE_DEST_FILES_LEN}-1")

    foreach(FILE_IDX RANGE ${FILES_IDX_END})

      list(GET PARSE_SOURCE_FILES ${FILE_IDX} SOURCE_FILE)
      set(SOURCE_FILE_FULL "${PARSE_SOURCE_DIR}/${PARSE_SOURCE_PREFIX}${SOURCE_FILE}")

      list(GET PARSE_DEST_FILES ${FILE_IDX} DEST_FILE)
      set(DEST_FILE_FULL "${PARSE_DEST_DIR}/${DEST_FILE}")

      #print_var(SOURCE_FILE_FULL)
      #print_var(DEST_FILE_FULL)

      add_custom_command(
        OUTPUT ${DEST_FILE_FULL}
        DEPENDS ${SOURCE_FILE_FULL}
        COMMAND ${CMAKE_COMMAND} ARGS -E copy ${SOURCE_FILE_FULL} ${DEST_FILE_FULL}
        )

      list(APPEND DEST_FILES_LIST "${DEST_FILE_FULL}")

    endforeach()
  else()
    foreach(SOURCE_FILE ${PARSE_SOURCE_FILES})
      set(DEST_FILE "${SOURCE_FILE}")

      set(SOURCE_FILE_FULL "${PARSE_SOURCE_DIR}/${PARSE_SOURCE_PREFIX}${SOURCE_FILE}")
      set(DEST_FILE_FULL "${PARSE_DEST_DIR}/${DEST_FILE}")

      #print_var(SOURCE_FILE_FULL)
      #print_var(DEST_FILE_FULL)

      add_custom_command(
        OUTPUT ${DEST_FILE_FULL}
        DEPENDS ${SOURCE_FILE_FULL}
        COMMAND ${CMAKE_COMMAND} ARGS -E copy ${SOURCE_FILE_FULL} ${DEST_FILE_FULL}
        )

      list(APPEND DEST_FILES_LIST "${DEST_FILE_FULL}")

    endforeach()
  endif()

  #print_var(DEST_FILES_LIST)

  if (PACKAGE_NAME AND NOT PARSE_NOEXEPREFIX)
    set(PACKAGE_PREFIX "${PACKAGE_NAME}_")
    set(TARGET_NAME "${PACKAGE_PREFIX}${TARGET_NAME}")
  endif()

  add_custom_target( ${TARGET_NAME} ALL
    DEPENDS ${DEST_FILES_LIST}
    )

  if (PARSE_EXEDEPS)
    foreach(EXEDEP ${PARSE_EXEDEPS})
      add_dependencies(${PACKAGE_PREFIX}${EXEDEP} ${TARGET_NAME})
    endforeach()
  endif()

  if (PARSE_TARGETDEPS)
    foreach(TARGETDEP ${PARSE_TARGETDEPS})
      add_dependencies(${TARGETDEP} ${TARGET_NAME})
    endforeach()
  endif()

endfunction()
