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


INCLUDE(TribitsAddTestHelpers)
INCLUDE(ParseVariableArguments)


#
# Function that copies a list of source files to a list of destination
# files between two base directories.
#
# This function has a few valid calling modes:
#
# 1) Source files and destination files have the same name:
#
#   TRIBITS_COPY_FILES_TO_BINARY_DIR(<TARGET_NAME>
#     SOURCE_FILES file1 file2 ...
#     [SOURCE_DIR sourceDir]
#     [DEST_DIR destDir]
#     [TARGETDEPS targDep1 targDep2 ...]
#     [EXEDEPS exeDep1 exeDep2 ...]
#     [NOEXEPREFIX]
#     [CATEGORIES <category1>  <category2> ...]
#     )
#
#   In this case, the names of the source files and the destination files
#   are the same but just live in different directories.
#
# 2) Source files have a prefix different from the destination files:
#
#   TRIBITS_COPY_FILES_TO_BINARY_DIR(<TARGET_NAME>
#     DEST_FILES file1 file2 ...
#     SOURCE_PREFIX srcPrefix
#     [SOURCE_DIR sourceDir]
#     [DEST_DIR destDir]
#     [EXEDEPS exeDep1 exeDep2 ...]
#     [NOEXEPREFIX]
#     [CATEGORIES <category1>  <category2> ...]
#     )
#
#   In this case, the source files have the same basic name as the
#   destination files except they have the prefix 'srcPrefix' appended
#   to the name.
#
# 3) Source files and destination files have completely different names:
#
#   TRIBITS_COPY_FILES_TO_BINARY_DIR(<TARGET_NAME>
#     SOURCE_FILES sfile1 sfile2 ...
#     [SOURCE_DIR sourceDir]
#     DEST_FILES dfile1 dfile2 ...
#     [DEST_DIR destDir]
#     [EXEDEPS exeDep1 exeDep2 ...]
#     [NOEXEPREFIX]
#     [CATEGORIES <category1>  <category2> ...]
#     )
#
#   In this case, the source files and destination files have completely
#   different prefixes.
#
# Here is the documentation for the individual arguments:
#
#   SOURCE_FILES file1 file2 ...
#     Listing of the source files relative to the source directory given
#     by the SOURCE_DIR argument 'sourceDir'.  If omited, this list will
#     be the same as DEST_FILES with the SOURCE_PREFIX argument 'srcPrefix'
#     appended.
#
#   SOURCE_DIR sourceDir
#     Optional argument that gives (absolute) the base directory for all of the
#     source files.  If omited, this takes the default value of 
#     ${CMAKE_CURRENT_SOURCE_DIR}.
#
#   DEST_FILES file1 file2 ...
#     Listing of the destination files relative to the destination directory given
#     by the DEST_DIR argument 'destDir'  If omited, this list will be
#     the same as given by the SOURCE_FILES list.
#
#   DEST_DIR destDir
#     Optional argument that gives the (absolute) base directory for all of the
#     destination files.  If omited, this takes the default value of 
#     ${CMAKE_CURRENT_BINARY_DIR}
#
#   TARGETDEPS targDep1 targDep2 ...
#     Listing of general CMake targets that these files will be added as
#     dependencies to.
#
#   EXEDEPS exeDep1 exeDep2 ...
#     Listing of executable targets that these files will be added as
#     dependencies to.  By default the prefix '${PACKAGE_NAME}_' will
#     is appended to the names of the targets.  Also, the '.exe' suffix
#     will be post-pended as well.  This ensures that if the executable
#     target is built that these files will also be copied as well.
#
#   NOEXEPREFIX
#     Option that determines if the prefix '${PACKAGE_NAME}_' will be appended
#     to the arguments in the EXEDEPS list. 
#

FUNCTION(TRIBITS_COPY_FILES_TO_BINARY_DIR TARGET_NAME)

  #
  # A) Parse input arguments
  #

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    "SOURCE_DIR;SOURCE_FILES;SOURCE_PREFIX;DEST_DIR;DEST_FILES;EXEDEPS;TARGETDEPS;CATEGORIES"
    #options
    "NOEXEPREFIX"
    ${ARGN}
    )

  SET(ADD_THE_TEST FALSE)
  TRIBITS_ADD_TEST_PROCESS_CATEGORIES(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  IF (NOT PARSE_SOURCE_DIR)
    SET(PARSE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  ENDIF()

  IF (NOT PARSE_DEST_DIR)
    SET(PARSE_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()

  SET(RENAME 1)

  IF (PARSE_SOURCE_FILES AND NOT PARSE_DEST_FILES)
    SET(PARSE_DEST_FILES ${PARSE_SOURCE_FILES})
    SET(RENAME 0)
  ENDIF()

  IF (PARSE_DEST_FILES AND NOT PARSE_SOURCE_FILES)
    SET(PARSE_SOURCE_FILES ${PARSE_DEST_FILES})
    SET(RENAME 0)
  ENDIF()

  #
  # B) Validate arguments
  #

  LIST(LENGTH PARSE_SOURCE_DIR PARSE_SOURCE_DIR_LEN)
  LIST(LENGTH PARSE_SOURCE_FILES PARSE_SOURCE_FILES_LEN)
  LIST(LENGTH PARSE_SOURCE_PREFIX PARSE_SOURCE_PREFIX_LEN)
  LIST(LENGTH PARSE_DEST_DIR PARSE_DEST_DIR_LEN)
  LIST(LENGTH PARSE_DEST_FILES PARSE_DEST_FILES_LEN)

  IF (PARSE_SOURCE_DIR_LEN GREATER 1)
    MESSAGE(SEND_ERROR "Error, there can only be 0 or one SOURCE_DIR arguments!")
  ENDIF()

  IF (PARSE_DEST_DIR_LEN GREATER 1)
    MESSAGE(SEND_ERROR "Error, there can only be 0 or one DEST_DIR arguments!")
  ENDIF()

  IF (PARSE_SOURCE_PREFIX_LEN GREATER 1)
    MESSAGE(SEND_ERROR "Error, If SOURCE_PREFIX can only take one argument!")
  ENDIF()

  IF (PARSE_SOURCE_FILES_LEN EQUAL 0)
    MESSAGE(SEND_ERROR "Error, there are no source files listed!")
  ENDIF()

  IF (PARSE_DEST_FILES_LEN EQUAL 0)
    MESSAGE(SEND_ERROR "Error, there are no destination files listed!")
  ENDIF()

  IF (NOT PARSE_SOURCE_FILES_LEN EQUAL ${PARSE_DEST_FILES_LEN})
    MESSAGE(SEND_ERROR "Error, there are not the same number of source files ${PARSE_SOURCE_FILES_LEN} and dest files ${PARSE_DEST_FILES_LEN}!")
  ENDIF()

  #
  # C) Build the list of command and dependencies
  #

  SET(DEST_FILES_LIST)
  IF(RENAME)
    MATH(EXPR FILES_IDX_END "${PARSE_DEST_FILES_LEN}-1")

    FOREACH(FILE_IDX RANGE ${FILES_IDX_END})

      LIST(GET PARSE_SOURCE_FILES ${FILE_IDX} SOURCE_FILE)
      SET(SOURCE_FILE_FULL "${PARSE_SOURCE_DIR}/${PARSE_SOURCE_PREFIX}${SOURCE_FILE}")

      LIST(GET PARSE_DEST_FILES ${FILE_IDX} DEST_FILE)
      SET(DEST_FILE_FULL "${PARSE_DEST_DIR}/${DEST_FILE}")

      #PRINT_VAR(SOURCE_FILE_FULL)
      #PRINT_VAR(DEST_FILE_FULL)

      ADD_CUSTOM_COMMAND(
        OUTPUT ${DEST_FILE_FULL}
        DEPENDS ${SOURCE_FILE_FULL}
        COMMAND ${CMAKE_COMMAND} ARGS -E copy ${SOURCE_FILE_FULL} ${DEST_FILE_FULL}
        )

      LIST(APPEND DEST_FILES_LIST "${DEST_FILE_FULL}")

    ENDFOREACH()
  ELSE()
    FOREACH(SOURCE_FILE ${PARSE_SOURCE_FILES})
      SET(DEST_FILE "${SOURCE_FILE}")

      SET(SOURCE_FILE_FULL "${PARSE_SOURCE_DIR}/${PARSE_SOURCE_PREFIX}${SOURCE_FILE}")
      SET(DEST_FILE_FULL "${PARSE_DEST_DIR}/${DEST_FILE}")

      #PRINT_VAR(SOURCE_FILE_FULL)
      #PRINT_VAR(DEST_FILE_FULL)

      ADD_CUSTOM_COMMAND(
        OUTPUT ${DEST_FILE_FULL}
        DEPENDS ${SOURCE_FILE_FULL}
        COMMAND ${CMAKE_COMMAND} ARGS -E copy ${SOURCE_FILE_FULL} ${DEST_FILE_FULL}
        )

      LIST(APPEND DEST_FILES_LIST "${DEST_FILE_FULL}")

    ENDFOREACH()
  ENDIF()

  #PRINT_VAR(DEST_FILES_LIST)

  IF (PACKAGE_NAME AND NOT PARSE_NOEXEPREFIX)
    SET(PACKAGE_PREFIX "${PACKAGE_NAME}_")
    SET(TARGET_NAME "${PACKAGE_PREFIX}${TARGET_NAME}")    
  ENDIF()

  ADD_CUSTOM_TARGET( ${TARGET_NAME} ALL
    DEPENDS ${DEST_FILES_LIST}
    )

  IF (PARSE_EXEDEPS)
    FOREACH(EXEDEP ${PARSE_EXEDEPS})
      ADD_DEPENDENCIES(${PACKAGE_PREFIX}${EXEDEP} ${TARGET_NAME})
    ENDFOREACH()
  ENDIF()

  IF (PARSE_TARGETDEPS)
    FOREACH(TARGETDEP ${PARSE_TARGETDEPS})
      ADD_DEPENDENCIES(${TARGETDEP} ${TARGET_NAME})
    ENDFOREACH()
  ENDIF()

ENDFUNCTION()
