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


#
# @MACRO: TRIBITS_INCLUDE_DIRECTORIES()
#
# This function is to override the standard behavior of the built-in CMake
# ``INCLUDE_DIRECTORIES()`` command.
#
# Usage::
#
#   TRIBITS_INCLUDE_DIRECTORIES(
#     [REQUIRED_DURING_INSTALLATION_TESTING] <dir0> <dir1> ...
#     )
#
# If specified, ``REQUIRED_DURING_INSTALLATION_TESTING`` can appear anywhere
# in the argument list.
#
# This function allows overriding the default behavior of
# ``INCLUDE_DIRECTORIES()`` for installation testing, to ensure that include
# directories will not be inadvertently added to the build lines for tests
# during installation testing (see `Installation and Backward Compatibility
# Testing`_). Normally we want the include directories to be handled as cmake
# usually does.  However during TriBITS installation testing we do not want
# most of the include directories to be used as the majority of the files
# should come from the installation we are building against.  There is an
# exception to this and that is when there are test only headers that are
# needed.  For that case ``REQUIRED_DURING_INSTALLATION_TESTING`` must be
# passed in to ensure the include paths are added for installation testing.
#
MACRO(TRIBITS_INCLUDE_DIRECTORIES)

  CMAKE_PARSE_ARGUMENTS(
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

  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_REQUIRED_DURING_INSTALLATION_TESTING)
    _INCLUDE_DIRECTORIES(${PARSE_UNPARSED_ARGUMENTS})
  ENDIF()
ENDMACRO()


# Deprecated.  Use TRIBITS_INCLUDE_DIRECTORIES() instead!
MACRO(INCLUDE_DIRECTORIES)

 CMAKE_PARSE_ARGUMENTS(
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

#  IF (PARSE_REQUIRED_DURING_INSTALLATION_TESTING)
#    MESSAGE(WARNING "Warning: the override INCLUDE_DIRECTORIES() is deprecated,"
#    " use TRIBITS_INCLUDE_DIRECTORIES() instead!")
#  ENDIF()

  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_REQUIRED_DURING_INSTALLATION_TESTING)
    _INCLUDE_DIRECTORIES(${PARSE_UNPARSED_ARGUMENTS})
  ENDIF()
ENDMACRO()
