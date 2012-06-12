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


INCLUDE(TribitsConstants)

INCLUDE(PrintVar)
INCLUDE(Split)


#
# Macro that processes the list of TPLs
#
# This macro reads from the varible
# ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS for a given repository and
# and fills the variables ${PROJECT_NAME}_TPLS, ${PROJECT_NAME}_NUM_TPLS,
# ${PROJECT_NAME}_REVERSE_TPLS.  For each TPL, it also sets the varible
# ${TPL_NAME}_FINDMOD and ${TPL_NAME}_CLASSIFICATION.
#

MACRO(TRIBITS_PROCESS_TPLS_LISTS  REPOSITORY_NAME  REPOSITORY_DIR)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("TRIBITS_PROCESS_TPLS_LISTS:  '${REPOSITORY_NAME}'  '${REPOSITORY_DIR}'")
  ENDIF()

  #SET(TRIBITS_PROCESS_TPLS_LISTS_DEBUG ON)
  SET(TRIBITS_PROCESS_TPLS_LISTS_DEBUG OFF)

  SET(TPL_NAME_OFFSET 0)
  SET(TPL_FINDMOD_OFFSET 1)
  SET(TPL_CLASSIFICATION_OFFSET 2)
  SET(TPL_NUM_COLUMNS 3)

  LIST(LENGTH ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS
    ${REPOSITORY_NAME}_CURR_NUM_TPLS_FULL)
  MATH(EXPR ${REPOSITORY_NAME}_CURR_NUM_TPLS
    "${${REPOSITORY_NAME}_CURR_NUM_TPLS_FULL}/${TPL_NUM_COLUMNS}")

  IF (${REPOSITORY_NAME}_CURR_NUM_TPLS GREATER 0)

    MATH(EXPR ${REPOSITORY_NAME}_LAST_TPL_IDX
      "${${REPOSITORY_NAME}_CURR_NUM_TPLS}-1")
  
    FOREACH(TPL_IDX RANGE ${${REPOSITORY_NAME}_LAST_TPL_IDX})

      IF (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        PRINT_VAR(TPL_IDX)
      ENDIF()
 
      # Get fields for this TPL
  
      MATH(EXPR TPL_NAME_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_NAME_OFFSET}")
      LIST(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_NAME_IDX}
        TPL_NAME)
      IF (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        PRINT_VAR(TPL_NAME)
      ENDIF()
  
      MATH(EXPR TPL_FINDMOD_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_FINDMOD_OFFSET}")
      LIST(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_FINDMOD_IDX}
        TPL_FINDMOD)
      IF (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        PRINT_VAR(TPL_FINDMOD)
      ENDIF()
  
      MATH(EXPR TPL_CLASSIFICATION_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_CLASSIFICATION_OFFSET}")
      LIST(GET ${REPOSITORY_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_CLASSIFICATION_IDX}
        TPL_CLASSIFICATION)
      IF (TRIBITS_PROCESS_TPLS_LISTS_DEBUG)
        PRINT_VAR(TPL_CLASSIFICATION)
      ENDIF()
  
      # Update TPLS list (unless the TPL already exists)
   
      IF (${TPL_NAME}_FINDMOD)
        # If the varaible ${TPL_NAME}_FINDMOD already exists, then this TPL
        # has already been defined in a previous repository.  In this case, we
        # will just leave the TPL in its current position.
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE("-- " "NOTE: The TPL ${TPL_NAME} has already been defined so leaving it"
            " in the same location and not adding it again!") 
        ENDIF()
      ELSE()
        LIST(APPEND ${PROJECT_NAME}_TPLS ${TPL_NAME})
      ENDIF() 
 
      # Set ${TPL_NAME}_CLASSIFICATION
  
      IF (TPL_CLASSIFICATION STREQUAL PS
        OR TPL_CLASSIFICATION STREQUAL SS
        OR TPL_CLASSIFICATION STREQUAL TS
        OR TPL_CLASSIFICATION STREQUAL EX
        )
      ELSE()
        MESSAGE(FATAL_ERROR "Error the TPL classification '${TPL_CLASSIFICATION}'"
          " for the TPL ${TPL_NAME} is not a valid classification." )
      ENDIF()
  
      IF (NOT ${TPL_NAME}_CLASSIFICATION) # Allow for testing override
        SET(${TPL_NAME}_CLASSIFICATION ${TPL_CLASSIFICATION})
      ENDIF()
  
      # Set ${TPL_NAME}_FINDMOD

      #PRINT_VAR(REPOSITORY_DIR)

      IF ("${REPOSITORY_DIR}" STREQUAL "." OR IS_ABSOLUTE ${TPL_FINDMOD})
        SET(REPOSITORY_DIR_AND_SEP "")
      ELSE()
        SET(REPOSITORY_DIR_AND_SEP "${REPOSITORY_DIR}/")
      ENDIF()
      #PRINT_VAR(REPOSITORY_DIR_AND_SEP)

      SET(TPL_FINDMOD "${REPOSITORY_DIR_AND_SEP}${TPL_FINDMOD}")
      #PRINT_VAR(TPL_FINDMOD)

      SET(TPL_FINDMOD_STD_NAME "FindTPL${TPL_NAME}.cmake")

      IF (TPL_FINDMOD)
        STRING(REGEX MATCH ".+/$" FINDMOD_IS_DIR "${TPL_FINDMOD}")
        #PRINT_VAR(FINDMOD_IS_DIR)
        IF (FINDMOD_IS_DIR)
          SET(${TPL_NAME}_FINDMOD "${TPL_FINDMOD}${TPL_FINDMOD_STD_NAME}")
        ELSE()
          SET(${TPL_NAME}_FINDMOD ${TPL_FINDMOD})
        ENDIF()
      ELSE()
        SET(${TPL_NAME}_FINDMOD ${TPL_FINDMOD_STD_NAME})
      ENDIF()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(${TPL_NAME}_FINDMOD)
      ENDIF()
  
      # Set the enable cache variable for ${TPL_NAME}
  
      MULTILINE_SET(DOCSTR
        "Enable support for the TPL ${TPL_NAME} in all supported ${PROJECT_NAME} packages."
        "  This can be set to 'ON', 'OFF', or left empty ''."
        )
      SET_CACHE_ON_OFF_EMPTY( TPL_ENABLE_${TPL_NAME} "" ${DOCSTR} )
  
      # 2008/11/25: rabartl: Above, we use the prefix TPL_ instead of
      # ${PROJECT_NAME}_ in order to make it clear that external TPLs are
      # different from packages so users don't get confused and
      # think that the project actually includes some TPL when it does not!
  
    ENDFOREACH()

  ENDIF()
  
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PROJECT_NAME}_TPLS)
  ENDIF()

  # Get the final length

  LIST(LENGTH ${PROJECT_NAME}_TPLS ${PROJECT_NAME}_NUM_TPLS)
  PRINT_VAR(${PROJECT_NAME}_NUM_TPLS)
  
  # Create a reverse list for later use
  
  IF (${PROJECT_NAME}_TPLS)
    SET(${PROJECT_NAME}_REVERSE_TPLS ${${PROJECT_NAME}_TPLS})
    LIST(REVERSE ${PROJECT_NAME}_REVERSE_TPLS)
  ELSE()
    SET(${PROJECT_NAME}_REVERSE_TPLS)
  ENDIF()

ENDMACRO()
