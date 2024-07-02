# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(CMakeParseArguments)
include(TribitsDeprecatedHelpers)


# @MACRO: tribits_include_directories()
#
# This function overrides the standard behavior of the built-in CMake
# ``include_directories()`` command for special behavior for installation
# testing.
#
# Usage::
#
#   tribits_include_directories(
#     [REQUIRED_DURING_INSTALLATION_TESTING] <dir0> <dir1> ... )
#
# If specified, ``REQUIRED_DURING_INSTALLATION_TESTING`` can appear anywhere
# in the argument list.
#
# This function allows overriding the default behavior of
# ``include_directories()`` for installation testing, to ensure that include
# directories will not be inadvertently added to the build lines for tests
# during installation testing (see `Installation and Backward Compatibility
# Testing`_). Normally we want the include directories to be handled as cmake
# usually does.  However during TriBITS installation testing we do not want
# most of the include directories to be used as the majority of the files
# should come from the installation we are building against.  The exception is
# when there are test only headers that are needed.  For that case
# ``REQUIRED_DURING_INSTALLATION_TESTING`` must be passed in to ensure the
# include paths are added for installation testing.
#
macro(tribits_include_directories)

  cmake_parse_arguments(
    #prefix
    PARSE
    #options
    "REQUIRED_DURING_INSTALLATION_TESTING"
    #one_value_keywords
    ""
    #mulit_value_keywords
    ""
    ${ARGN}
    )

  if(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
      OR PARSE_REQUIRED_DURING_INSTALLATION_TESTING
    )
    if (TRIBITS_HIDE_DEPRECATED_INCLUDE_DIRECTORIES_OVERRIDE)
      include_directories(${PARSE_UNPARSED_ARGUMENTS})
    else()
      _include_directories(${PARSE_UNPARSED_ARGUMENTS})
    endif()
  endif()
endmacro()


if (NOT TRIBITS_HIDE_DEPRECATED_INCLUDE_DIRECTORIES_OVERRIDE)

# Deprecated.  Use tribits_include_directories() instead!
#
# To hide this macro from even being defined, set
# ``TRIBITS_HIDE_DEPRECATED_INCLUDE_DIRECTORIES_OVERRIDE=TRUE``.
#
macro(include_directories)

  tribits_deprecated_command(include_directories
    MESSAGE "Use tribits_include_directories() instead."
    )

  cmake_parse_arguments(
    #prefix
    PARSE
    #options
    "REQUIRED_DURING_INSTALLATION_TESTING"
    #one_value_keywords
    ""
    #mulit_value_keywords
    ""
    ${ARGN}
    )

  if(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
      OR PARSE_REQUIRED_DURING_INSTALLATION_TESTING
    )
    _include_directories(${PARSE_UNPARSED_ARGUMENTS})
  endif()
endmacro()

endif (NOT TRIBITS_HIDE_DEPRECATED_INCLUDE_DIRECTORIES_OVERRIDE)
