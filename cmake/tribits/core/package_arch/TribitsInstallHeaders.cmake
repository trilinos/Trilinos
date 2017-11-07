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
# @FUNCTION: TRIBITS_INSTALL_HEADERS()
#
# Function used to (optionally) install header files using ``INSTALL()``
# command.
#
# Usage::
#
#   TRIBITS_INSTALL_HEADERS(
#     HEADERS <h0> <h1> ...
#     [INSTALL_SUBDIR <subdir>]
#     [COMPONENT <component>]
#     )
#
# The formal arguments are:
#
#   ``HEADERS <h0> <h1> ...``
#
#     List of header files to install.  By default, these header files are
#     assumed to be in the current source directory.  They can also contain
#     the relative path or absolute path to the files if they are not in the
#     current source directory.
#
#   ``INSTALL_SUBDIR <subdir>``
#
#     Optional subdirectory that the headers will be installed under the
#     standard installation directory.  If ``<subdir>!=""``, then the headers
#     will be installed under
#     ``${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/<subdir>``.  Otherwise, they will
#     be installed under ``${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/``.
#
#   ``COMPONENT <component>``
#
#     If specified, then ``COMPONENT <component>`` will be passed into
#     ``INSTALL()``.  Otherwise, ``COMPONENT ${PROJECT_NAME}`` will get used.
#
# If `${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS`_ is ``FASLE``, then the
# headers will not get installed.
#
FUNCTION(TRIBITS_INSTALL_HEADERS)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    SET(TRIBITS_INSTALL_HEADERS_DEBUG_DUMP  TRUE)
  ENDIF()

  IF (TRIBITS_INSTALL_HEADERS_DEBUG_DUMP)
    MESSAGE("\nTRIBITS_INSTALL_HEADERS: ${ARGN}")
  ENDIF()

  CMAKE_PARSE_ARGUMENTS(
    #prefix
    PARSE
    #Options
    ""
    #one_value_keywords
    ""
    #multi_value_keywords
    "HEADERS;INSTALL_SUBDIR;COMPONENT"
    ${ARGN}
    )

  TRIBITS_CHECK_FOR_UNPARSED_ARGUMENTS()

  # ToDo: Assert PARSE_HEADERS has at least one argument!
  # ToDo: Assert PARSE_INSTALL_DIR has 0 or 1 argumnets!
  # ToDo: Assert PARSE_COMONENT has 0 or 1 argumnets!
  
  IF (PARSE_INSTALL_SUBDIR)
    SET(INSTALL_DIR "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/${PARSE_INSTALL_SUBDIR}")
  ELSE()
    SET(INSTALL_DIR "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")
  ENDIF()

  IF (PARSE_COMPONENT)
    SET(COMPONENT ${PARSE_COMPONENT})
  ELSE()
    SET(COMPONENT ${PROJECT_NAME})
  ENDIF()

  IF (${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS)
    IF (TRIBITS_INSTALL_HEADERS_DEBUG_DUMP)
      MESSAGE("\nTRIBITS_INSTALL_HEADERS: Installing headers into '${INSTALL_DIR}'")
    ENDIF()
    INSTALL(
      FILES ${PARSE_HEADERS}
      DESTINATION "${INSTALL_DIR}"
      COMPONENT ${COMPONENT}
      )
  ENDIF()  

ENDFUNCTION()
