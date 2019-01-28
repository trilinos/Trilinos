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


INCLUDE(TribitsHostType)


#
# This module defines the datastructure for the list of packages
# ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS which has the form:
#
#   Package0_name  Package0_dir Package0_classification
#   Package1_name  Package1_dir Package1_classification
#   ...
#
# There are 3 fields per row all stored in a flat array.
#


SET(PLH_NUM_FIELDS_PER_PACKAGE 3)
SET(PLH_NUM_PACKAGE_DIR_OFFSET 1)
SET(PLH_NUM_PACKAGE_CLASSIFICATION_OFFSET 2)


MACRO(TRIBITS_SET_PACKAGE_TO_EX  PACKAGE_NAME)
  LIST(FIND ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
    ${PACKAGE_NAME} PACKAGE_NAME_IDX)
  IF (PACKAGE_NAME_IDX EQUAL -1)
    MESSAGE(
      "\n***"
      "\n*** NOTE: Package ${PACKAGE_NAME} not found in list of packages!"
      "\n***\n"
      )
  ELSE()
    MATH(EXPR PACKAGE_CLASSIFICATION_IDX
      "${PACKAGE_NAME_IDX}+${PLH_NUM_PACKAGE_CLASSIFICATION_OFFSET}")
    LIST(INSERT ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
      ${PACKAGE_CLASSIFICATION_IDX} EX )
    MATH(EXPR PACKAGE_CLASSIFICATION_IDX_PLUS_1
      "${PACKAGE_CLASSIFICATION_IDX}+1")
    LIST(REMOVE_AT ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
      ${PACKAGE_CLASSIFICATION_IDX_PLUS_1} )
  ENDIF()

ENDMACRO()


#
# @MACRO: TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS()
#
# Disable a package automatically for a list of platforms.
#
# Usage::
#
#   TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS( <packageName>
#     <hosttype0> <hosttype1> ...)
#
# If any of the host-type arguments ``<hosttypei>`` matches the
# ``${PROJECT_NAME}_HOSTTYPE`` variable for the current platform, then package
# ``<packageName>`` test group classification is changed to ``EX``.  Changing
# the package test group classification to ``EX`` results in the package being
# disabled by default (see `EX SE packages disabled by default`_).  However,
# an explicit enable can still enable the package.
#
MACRO( TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS  PACKAGE_NAME )
  #MESSAGE("TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS: ${PACKAGE_NAME}")
  #PRINT_VAR(${PROJECT_NAME}_HOSTTYPE)
  FOREACH(HOSTTYPE ${ARGN})
    #PRINT_VAR(HOSTTYPE)
    IF (${PROJECT_NAME}_HOSTTYPE STREQUAL ${HOSTTYPE})
      #MESSAGE("${${PROJECT_NAME}_HOSTTYPE} == ${HOSTTYPE}")
      TRIBITS_SET_PACKAGE_TO_EX(${PACKAGE_NAME})
      #PRINT_VAR(${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)
      IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
        MESSAGE(
          "\n***"
          "\n*** NOTE: User has set ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON but the"
          "\n*** package ${PACKAGE_NAME} is not supported on this platform type '${HOSTTYPE}'!"
          "\n***\n"
          )
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDMACRO()


MACRO( PACKAGE_DISABLE_ON_PLATFORMS  PACKAGE_NAME_IN_ )
  MESSAGE(WARNING "PACKAGE_DISABLE_ON_PLATFORMS() is deprecated!"
    "  Use TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS() instead!")
  TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(${PACKAGE_NAME_IN_} ${ARGN})
ENDMACRO()



FUNCTION(TRIBITS_UPDATE_PS_PT_SS_ST  THING_TYPE  THING_NAME  TESTGROUP_VAR)

  #MESSAGE("TRIBITS_UPDATE_PS_PT_SS_ST:  ${THING_TYPE}  ${THING_NAME}  ${TESTGROUP_VAR}")

  SET(TESTGROUP_IN ${${TESTGROUP_VAR}})
  #PRINT_VAR(TESTGROUP_IN)

  IF (TESTGROUP_IN STREQUAL PS)
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "WARNING: ${THING_TYPE} ${THING_NAME} TESTGROUP 'PS' is depricated."
        "  Use 'PT' instead!")
    ENDIF()
    SET(TESTGROUP_OUT PT)
  ELSEIF (TESTGROUP_IN STREQUAL SS)
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "WARNING: ${THING_TYPE} ${THING_NAME} TESTGROUP 'SS' is depricated."
        "  Use 'ST' instead!")
    ENDIF()
    SET(TESTGROUP_OUT ST)
  ELSEIF (TESTGROUP_IN STREQUAL TS)
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "WARNING: ${THING_TYPE} ${THING_NAME} TESTGROUP 'TS' is depricated."
        "  Use 'TT' instead!")
    ENDIF()
    SET(TESTGROUP_OUT TT)
  ELSE()
    SET(TESTGROUP_OUT ${TESTGROUP_IN})
  ENDIF()
  #PRINT_VAR(TESTGROUP_OUT)

  SET(${TESTGROUP_VAR} ${TESTGROUP_OUT} PARENT_SCOPE)

ENDFUNCTION()
