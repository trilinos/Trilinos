# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(CMakeParseArguments)
include(PrintVar)


#
# Function that creates a set of template code client headers that either
# includes or does not include the template instantiations depending on
# whether implicit or explicit instantiation is supported or not.
#
# Usage:
#
#    tribits_create_client_template_headers(
#      BASE_DIR
#      [ADDITIONAL_OUTPUT_DIRS ABSDIR1 ABSDIR2 ...]
#      )
#
# The arguments are:
#
#  BASE_DIR
#
#    The base directory where files with the extension
#    ${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT} will be
#    globed for.
#
#  ADDITIONAL_OUTPUT_DIRS
#
#    If set, then the files will be copied to an additional output
#    directories as well.  These must be absolute paths.
#
# The default file extensions are:
#
#    ${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT = "_decl.hpp"
#    ${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT = "_def.hpp"
#

function(tribits_create_client_template_headers BASE_DIR)

  #print_var(BASE_DIR)

  #
  # A) Parse the input arguments
  #

  cmake_parse_arguments(
    #prefix
    PARSE
    #options
    ""
    #one_value_keywords
    ""
    #multi_value_keywords
    "ADDITIONAL_OUTPUT_DIRS"
    ${ARGN}
    )

  tribits_check_for_unparsed_arguments()

  #
  # B) Get the names of the extensions
  #

  if (NOT ${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT)
    set(${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT "_decl.hpp")
  endif()

  if (NOT ${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT)
    set(${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT "_def.hpp")
  endif()

  #
  # C) Glob the names of all the X_decl.hpp files
  #

  file(GLOB DECL_HEADERS_LIST "${BASE_DIR}/*${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT}")
  #print_var(DECL_HEADERS_LIST)

  #
  # D) Write the client header files for each globed decl file
  #

  assert_defined(HAVE_${PARENT_PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION)

  foreach(DECL_HEADER ${DECL_HEADERS_LIST})

    # Get the base file names (without _decl.hpp)
    string(REGEX REPLACE ".*/(.+)${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT}" "\\1"  DECL_HEADER_BASE ${DECL_HEADER})
    #print_var(DECL_HEADER_BASE)

    # Create the client header file
    set(CLIENT_HEADER_STR "")
    string(APPEND CLIENT_HEADER_STR
      "#include \"${DECL_HEADER_BASE}${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT}\"\n"
       )
    if (HAVE_${PARENT_PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION)
        set(TEMPLATE_INSTANT_TYPE_NAME "explicit instantiation")
    else()
      set(TEMPLATE_INSTANT_TYPE_NAME "implicit instantiation")
       string(APPEND CLIENT_HEADER_STR
        "#include \"${DECL_HEADER_BASE}${${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT}\"\n"
         )
    endif()
    set(BIN_HEADER_FILE "${CMAKE_CURRENT_BINARY_DIR}/${DECL_HEADER_BASE}.hpp")
    set(WRITE_NEW_HEADER_FILE TRUE)
    if (EXISTS "${BIN_HEADER_FILE}")
      # See if the file is the same and if it is, skip writing it again to avoid
      # unnecessarily rebuilding object code.
      file(READ "${BIN_HEADER_FILE}" EXISTING_BIN_HEADER_STR)
      if (CLIENT_HEADER_STR STREQUAL EXISTING_BIN_HEADER_STR)
        set(WRITE_NEW_HEADER_FILE FALSE)
      endif()
    endif()

    # Write the header file
    if (WRITE_NEW_HEADER_FILE)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("Writing ${TEMPLATE_INSTANT_TYPE_NAME} header ${BIN_HEADER_FILE}")
      endif()
      file(WRITE "${BIN_HEADER_FILE}" "${CLIENT_HEADER_STR}")
    endif()

    # Create the SIERRA BJAM version of the header file
    foreach(OUTPUT_DIR ${PARSE_ADDITIONAL_OUTPUT_DIRS})
      set(EXTERNAL_CLIENT_HEADER_STR "")
      string(APPEND EXTERNAL_CLIENT_HEADER_STR
        "#include \"${DECL_HEADER_BASE}${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT}\"\n"
        "#ifndef HAVE_${PARENT_PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION\n"
        "#  include \"${DECL_HEADER_BASE}${${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT}\"\n"
        "#endif\n"
         )
      set(EXTERNAL_HEADER "${OUTPUT_DIR}/${DECL_HEADER_BASE}.hpp")
      if (NOT EXISTS "${EXTERNAL_HEADER}")
        file(WRITE "${EXTERNAL_HEADER}" "${EXTERNAL_CLIENT_HEADER_STR}")
      endif()
    endforeach()

  endforeach()

endfunction()

