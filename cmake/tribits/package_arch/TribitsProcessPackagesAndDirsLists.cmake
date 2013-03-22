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


INCLUDE(SetCacheOnOffEmpty)
INCLUDE(MultilineSet)
INCLUDE(AdvancedOption)
INCLUDE(AssertDefined)
INCLUDE(PrintVar)
INCLUDE(TribitsListHelpers)


#
# Helper functions
#


#
# Macro used in Dependencies.cmake files to allow some upstream dependent packages
# to be missing.
#
# Usage:
#
#   TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(<pack_1> <pack_2> ...)
#
# If the missing upstream SE package <pack_i> is optional, then the effect
# will be to simply ignore the missing package and remove it from the
# dependency list.  However, if the missing upstream SE package <pack_i> is
# required, then in addition to ignoring the missing package, the current SE
# (sub)package will also e hard disabled,
# i.e. ${PROJECT_NAME}_ENABLE_{CURRENT_PACKAGE}=OFF.
#
# This function is typically used in packages in external TriBITS repos that
#  depend on other packages in other exteral TriBITS repos that might be
#  missing.
#
# NOTE: Using this function effectively turns off error checking for
#  misspelled package names so it is important to only use it when it
#  absolutely is needed.
#
MACRO(TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES)
  FOREACH(TRIBITS_PACKAGE ${ARGN})
    SET(${TRIBITS_PACKAGE}_ALLOW_MISSING_EXTERNAL_PACKAGE TRUE)
  ENDFOREACH()
ENDMACRO()


#
# Below, we change the value of user cache values like
# ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME},
# ${PACKAGE_NAME}_ENABLE_TESTS, and ${PACKAGE_NAME}_ENABLE_EXAMPLES by
# just setting them to regular variables that live in the top scope
# instead of putting them back in the cache.  That means that they are
# used as global variables but we don't want to disturb the cache
# since that would change the behavior for future invocations of cmake
# (which is very confusing).  Because of this, these macros must all
# be called from the top-level ${PROJECT_NAME} CMakeLists.txt file and
# macros must call macros as not to change the variable scope.
#
# I had to do it this way in order to be able to get the right behavior which
# is:
#
# 1) Override the value of these variables in all CMake processing
#
# 2) Avoid changing the user cache values because that would be confusing and
# would make it hard to change package enables/disable later without blowing
# away the cache
# 


#
# Macro that sets up standard user options for each package
#

MACRO(TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS  PACKAGE_NAME  PACKAGE_CLASSIFICATION)

  #MESSAGE("TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS: ${PACKAGE_NAME} ${PACKAGE_CLASSIFICATION}")

  IF (${PACKAGE_CLASSIFICATION} STREQUAL PS OR ${PACKAGE_CLASSIFICATION} STREQUAL SS) 
    #MESSAGE("PS or SS")
    SET(PACKAGE_ENABLE "")
  ELSEIF (${PACKAGE_CLASSIFICATION} STREQUAL EX)
    #MESSAGE("EX")
    SET(PACKAGE_ENABLE OFF)
  ELSE()
    MESSAGE(FATAL_ERROR "Error the package classification '${PACKAGE_CLASSIFICATION}'"
      " for the package ${PACKAGE_NAME} is not a valid classification." )
  ENDIF()
  #PRINT_VAR(PACKAGE_ENABLE)

  IF (NOT ${PACKAGE_NAME}_CLASSIFICATION) # Allow testing override
    SET(${PACKAGE_NAME}_CLASSIFICATION "${PACKAGE_CLASSIFICATION}")
  ENDIF()

  MULTILINE_SET(DOCSTR
    "Enable the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave"
    " empty to allow for other logic to decide."
    )
  SET_CACHE_ON_OFF_EMPTY( ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
    "${PACKAGE_ENABLE}" ${DOCSTR} )

ENDMACRO()


#
# Function that determines if a package is to be excluded from implicit enable
#  as determined by it's project's
#  ${PARENT_REPO_NAME}_NO_IMPLICIT_PACKAGE_ENABLE[_EXCEPT] variables.
#

FUNCTION(TRIBITS_REPOSITORY_IMLICIT_PACKAGE_ENABLE_IS_ALLOWED
  UPSTREAM_PACKAGE_NAME_IN  PACKAGE_NAME_IN
  IMPLICIT_PACKAGE_ENABLE_ALLOWED_OUT
  )

  SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED TRUE)

  SET(PARENT_REPO_NAME ${${PACKAGE_NAME_IN}_PARENT_REPOSITORY})
  IF (IMPLICIT_PACKAGE_ENABLE_ALLOWED AND
    ${PARENT_REPO_NAME}_NO_IMPLICIT_PACKAGE_ENABLE
    )
    FIND_LIST_ELEMENT(${PARENT_REPO_NAME}_NO_IMPLICIT_PACKAGE_ENABLE_EXCEPT ${PACKAGE_NAME_IN}
      PACKAGE_EXCEPTED)
    IF (PACKAGE_EXCEPTED)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("-- " "NOTE: ${PACKAGE_NAME_IN} allowed to be implicitly enabled because "
          "${PARENT_REPO_NAME}_NO_IMPLICIT_PACKAGE_ENABLE=ON "
          "but package is included in the list ${PARENT_REPO_NAME}_NO_IMPLICIT_PACKAGE_ENABLE_EXCEPT!")
      ENDIF()
    ELSE()
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        IF (UPSTREAM_PACKAGE_NAME_IN)
          MESSAGE("-- " "WARNING: ${PACKAGE_NAME_IN} not being implicitly enabled even though"
            " ${UPSTREAM_PACKAGE_NAME_IN} has an optional dependence on it because"
            " ${PARENT_REPO_NAME}_NO_IMPLICIT_PACKAGE_ENABLE=ON!")
        ELSE()
          MESSAGE("-- " "WARNING: ${PACKAGE_NAME_IN} not being implicitly enabled because "
            "${PARENT_REPO_NAME}_NO_IMPLICIT_PACKAGE_ENABLE=ON!")
          ENDIF()
      ENDIF()
      SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED FALSE)
    ENDIF()
  ENDIF()

  SET(${IMPLICIT_PACKAGE_ENABLE_ALLOWED_OUT} ${IMPLICIT_PACKAGE_ENABLE_ALLOWED}
    PARENT_SCOPE )

ENDFUNCTION()


#
# Function that determines if it is okay to allow an implicit package enable
# based on its classification.
#

FUNCTION(TRIBITS_IMPLICIT_PACKAGE_ENABLE_IS_ALLOWED  UPSTREAM_PACKAGE_NAME_IN  PACKAGE_NAME_IN
  IMPLICIT_PACKAGE_ENABLE_ALLOWED_OUT
  )

  IF (${PACKAGE_NAME_IN}_CLASSIFICATION STREQUAL PS)
    SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED TRUE)
  ELSEIF (${PACKAGE_NAME_IN}_CLASSIFICATION STREQUAL SS
    AND ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE
    )
    SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED TRUE)
  ELSE()
    IF (UPSTREAM_PACKAGE_NAME_IN)
      MESSAGE("-- " "WARNING: Not Setting ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME_IN}=ON"
        " even though ${UPSTREAM_PACKAGE_NAME_IN} has an optional dependence on"
        " ${PACKAGE_NAME_IN} because ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE=OFF" )
    ENDIF()
    SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED FALSE)
  ENDIF()
 
  TRIBITS_REPOSITORY_IMLICIT_PACKAGE_ENABLE_IS_ALLOWED(
    "${UPSTREAM_PACKAGE_NAME_IN}"  "${PACKAGE_NAME_IN}"
    REPOSITORY_IMPLICIT_PACKAGE_ENABLE_ALLOWED_OUT
    )
  IF (NOT REPOSITORY_IMPLICIT_PACKAGE_ENABLE_ALLOWED_OUT)
    SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED FALSE)
  ENDIF()

  SET(${IMPLICIT_PACKAGE_ENABLE_ALLOWED_OUT} ${IMPLICIT_PACKAGE_ENABLE_ALLOWED}
    PARENT_SCOPE )

ENDFUNCTION()


#
# Macro that processes ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS into
# ${PROJECT_NAME}_PACKAGES, ${PROJECT_NAME}_PACKAGE_DIRS, ${PROJECT_NAME}_NUM_PACKAGES,
# ${PROJECT_NAME}_LAST_PACKAGE_IDX, and ${PROJECT_NAME}_REVERSE_PACKAGES.
#
# This macro also sets up the standard package options along with
# default enables/disables.
#
# NOTE: Set TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE=TRUE to see really verbose
# debug ouptut from this macro.
#

MACRO(TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS  REPOSITORY_NAME  REPOSITORY_DIR)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS:  '${REPOSITORY_NAME}'  '${REPOSITORY_DIR}'")
  ENDIF()

  #
  # Separate out separate lists of package names and directories
  #

  # Get the total number of packages defined  

  ASSERT_DEFINED(${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)
  IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
    PRINT_VAR(${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)
  ENDIF()
  LIST(LENGTH ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
    ${REPOSITORY_NAME}_NUM_PACKAGES_AND_FIELDS )
  MATH(EXPR ${REPOSITORY_NAME}_NUM_PACKAGES
    "${${REPOSITORY_NAME}_NUM_PACKAGES_AND_FIELDS}/${PLH_NUM_FIELDS_PER_PACKAGE}")
  IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
    PRINT_VAR(${REPOSITORY_NAME}_NUM_PACKAGES)
  ENDIF()
  MATH(EXPR ${REPOSITORY_NAME}_LAST_PACKAGE_IDX "${${REPOSITORY_NAME}_NUM_PACKAGES}-1")
 
  # Process each of the packages defined

  IF (${REPOSITORY_NAME}_NUM_PACKAGES GREATER 0)
  
    FOREACH(PACKAGE_IDX RANGE ${${REPOSITORY_NAME}_LAST_PACKAGE_IDX})
  
      IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        PRINT_VAR(${PROJECT_NAME}_PACKAGES)
      ENDIF()
  
      MATH(EXPR PACKAGE_NAME_IDX "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+0")
      LIST(GET ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
        ${PACKAGE_NAME_IDX} TRIBITS_PACKAGE )
      IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        PRINT_VAR(TRIBITS_PACKAGE)
      ENDIF()
  
      MATH(EXPR PACKAGE_DIR_IDX
        "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+${PLH_NUM_PACKAGE_DIR_OFFSET}")
      LIST(GET ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
        ${PACKAGE_DIR_IDX} PACKAGE_DIR )
      IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        PRINT_VAR(PACKAGE_DIR)
      ENDIF()
  
      MATH(EXPR PACKAGE_CLASSIFICATION_IDX
        "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+${PLH_NUM_PACKAGE_CLASSIFICATION_OFFSET}")
      LIST(GET ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
        ${PACKAGE_CLASSIFICATION_IDX} PACKAGE_CLASSIFICATION )
      IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        PRINT_VAR(PACKAGE_CLASSIFICATION)
      ENDIF()
  
      IF ("${REPOSITORY_DIR}" STREQUAL ".")
        SET(REPOSITORY_AND_PACKAGE_DIR "${PACKAGE_DIR}")
      ELSEIF("${PACKAGE_DIR}" STREQUAL ".")
        SET(REPOSITORY_AND_PACKAGE_DIR "${REPOSITORY_DIR}")
      ELSE()
        SET(REPOSITORY_AND_PACKAGE_DIR "${REPOSITORY_DIR}/${PACKAGE_DIR}")
      ENDIF()
      IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        PRINT_VAR(REPOSITORY_AND_PACKAGE_DIR)
      ENDIF()
  
      SET(PACKAGE_ABS_DIR "${PROJECT_SOURCE_DIR}/${REPOSITORY_AND_PACKAGE_DIR}")
      IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        PRINT_VAR(PROJECT_SOURCE_DIR)
        PRINT_VAR(REPOSITORY_AND_PACKAGE_DIR)
        PRINT_VAR(PACKAGE_ABS_DIR)
      ENDIF()
  
      IF (EXISTS ${PACKAGE_ABS_DIR})
        SET(PACKAGE_EXISTS TRUE)
      ELSE()
        SET(PACKAGE_EXISTS FALSE)
      ENDIF()
  
      IF (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES AND NOT PACKAGE_EXISTS)
        MESSAGE(
          "\n***"
          "\n*** Error, the package ${TRIBITS_PACKAGE} directory ${PACKAGE_ABS_DIR} does not exist!"
          "\n***\n" )
        MESSAGE(FATAL_ERROR "Stopping due to above error!")
      ENDIF()
  
      IF (PACKAGE_EXISTS OR ${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK)
        LIST(APPEND ${PROJECT_NAME}_PACKAGES ${TRIBITS_PACKAGE})
        LIST(APPEND ${PROJECT_NAME}_PACKAGE_DIRS "${REPOSITORY_AND_PACKAGE_DIR}")
        TRIBITS_INSERT_STANDARD_PACKAGE_OPTIONS(${TRIBITS_PACKAGE} ${PACKAGE_CLASSIFICATION})
        SET(${TRIBITS_PACKAGE}_PARENT_PACKAGE "")
        SET(${TRIBITS_PACKAGE}_PARENT_REPOSITORY ${REPOSITORY_NAME})
      ELSE()
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE(
            "\n***"
            "\n*** WARNING: Excluding package ${TRIBITS_PACKAGE} because ${PACKAGE_ABS_DIR}"
              " does not exist!"
            "\n***\n" )
        ENDIF()
      ENDIF()
  
      IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        PRINT_VAR(${TRIBITS_PACKAGE}_PARENT_PACKAGE)
        PRINT_VAR(${TRIBITS_PACKAGE}_PARENT_REPOSITORY)
      ENDIF()

      IF (TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS_VERBOSE)
        PRINT_VAR(${PROJECT_NAME}_PACKAGES)
      ENDIF()
  
    ENDFOREACH()
  
    # Get the actual number of packages that actually exist
  
    LIST(LENGTH ${PROJECT_NAME}_PACKAGES ${PROJECT_NAME}_NUM_PACKAGES )
    MATH(EXPR ${PROJECT_NAME}_LAST_PACKAGE_IDX "${${PROJECT_NAME}_NUM_PACKAGES}-1")
    
    # Create a reverse list for later use
    
    SET(${PROJECT_NAME}_REVERSE_PACKAGES ${${PROJECT_NAME}_PACKAGES})
    LIST(REVERSE ${PROJECT_NAME}_REVERSE_PACKAGES)

  ELSE()

    SET(${REPOSITORY_NAME}_NUM_PACKAGES 0)

  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${REPOSITORY_NAME}_NUM_PACKAGES)
  ENDIF()

  PRINT_VAR(${PROJECT_NAME}_NUM_PACKAGES)
  
  # Print the final set of packages in debug mode

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PROJECT_NAME}_PACKAGES)
    PRINT_VAR(${PROJECT_NAME}_PACKAGE_DIRS)
  ENDIF()

ENDMACRO()
