# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


INCLUDE(ParseVariableArguments)
INCLUDE(PrintVar)


#
# Function that creates a set of template code client headers that either
# includes or does not include the template instantiations depending on
# whether implicit or explicit instantiation is supported or not.
#
# Usage:
#
#    TRIBITS_CREATE_CLIENT_TEMPLATE_HEADERS(
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
#    directories as well.  These must be abolute paths.
#
# The default file extensions are:
#
#    ${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT = "_decl.hpp"
#    ${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT = "_def.hpp"
#

FUNCTION(TRIBITS_CREATE_CLIENT_TEMPLATE_HEADERS BASE_DIR)

  #PRINT_VAR(BASE_DIR)
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    "ADDITIONAL_OUPTUT_DIRS"
    #options
    ""
    ${ARGN}
    )

  #
  # B) Get the names of the extensions
  #

  IF (NOT ${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT)
    SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT "_decl.hpp")
  ENDIF()

  IF (NOT ${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT)
    SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT "_def.hpp")
  ENDIF()

  #
  # C) Glob the names of all the X_decl.hpp files
  #

  FILE(GLOB DECL_HEADERS_LIST "${BASE_DIR}/*${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT}")
  #PRINT_VAR(DECL_HEADERS_LIST)

  #
  # D) Write the client header files for each globed decl file
  #

  FOREACH(DECL_HEADER ${DECL_HEADERS_LIST})
 
    # Get the base file names (without _decl.hpp)
    STRING(REGEX REPLACE ".*/(.+)${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT}" "\\1"  DECL_HEADER_BASE ${DECL_HEADER})
    #PRINT_VAR(DECL_HEADER_BASE)

    # Create the client header file
    SET(CLIENT_HEADER_STR "")
    APPEND_STRING_VAR(CLIENT_HEADER_STR
      "#include \"${DECL_HEADER_BASE}${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT}\"\n"
       )
    IF (NOT HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION)
      APPEND_STRING_VAR(CLIENT_HEADER_STR
        "#include \"${DECL_HEADER_BASE}${${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT}\"\n"
         )
    ENDIF()
    #PRINT_VAR(CLIENT_HEADER_STR)
    SET(BIN_HEADER_FILE "${CMAKE_CURRENT_BINARY_DIR}/${DECL_HEADER_BASE}.hpp")
    SET(WRITE_NEW_HEADER_FILE TRUE)
    IF (EXISTS "${BIN_HEADER_FILE}")
      # See if the file is the same and if it is, skip writing it again to avoid
      # unecessarily rebuilding object code.
      FILE(READ "${BIN_HEADER_FILE}" EXISTING_BIN_HEADER_STR)
      IF (CLIENT_HEADER_STR STREQUAL EXISTING_BIN_HEADER_STR)
        SET(WRITE_NEW_HEADER_FILE FALSE)
      ENDIF()
    ENDIF()
    IF (WRITE_NEW_HEADER_FILE)
      FILE(WRITE "${BIN_HEADER_FILE}" "${CLIENT_HEADER_STR}")
    ENDIF()

    # Create the SIERRA BJAM version of the header file
    FOREACH(OUTPUT_DIR ${PARSE_ADDITIONAL_OUTPUT_DIRS})
      SET(EXTERNAL_CLIENT_HEADER_STR "")
      APPEND_STRING_VAR(EXTERNAL_CLIENT_HEADER_STR
        "#include \"${DECL_HEADER_BASE}${${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT}\"\n"
        "#ifndef HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION\n"
        "#  include \"${DECL_HEADER_BASE}${${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT}\"\n"
        "#endif\n"
         )
      SET(EXTERNAL_HEADER "${OUTPUT_DIR}/${DECL_HEADER_BASE}.hpp")
      IF (NOT EXISTS "${EXTERNAL_HEADER}")
        FILE(WRITE "${EXTERNAL_HEADER}" "${EXTERNAL_CLIENT_HEADER_STR}")
      ENDIF()
    ENDFOREACH()

  ENDFOREACH()

ENDFUNCTION()

