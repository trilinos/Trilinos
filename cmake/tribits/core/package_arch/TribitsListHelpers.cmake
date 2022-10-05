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


include(TribitsHostType)


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


set(PLH_NUM_FIELDS_PER_PACKAGE 3)
set(PLH_NUM_PACKAGE_DIR_OFFSET 1)
set(PLH_NUM_PACKAGE_CLASSIFICATION_OFFSET 2)


macro(tribits_set_package_to_ex  PACKAGE_NAME)
  list(FIND ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
    ${PACKAGE_NAME} PACKAGE_NAME_IDX)
  if (PACKAGE_NAME_IDX EQUAL -1)
    message(
      "\n***"
      "\n*** NOTE: Package ${PACKAGE_NAME} not found in list of packages!"
      "\n***\n"
      )
  else()
    math(EXPR PACKAGE_CLASSIFICATION_IDX
      "${PACKAGE_NAME_IDX}+${PLH_NUM_PACKAGE_CLASSIFICATION_OFFSET}")
    list(INSERT ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
      ${PACKAGE_CLASSIFICATION_IDX} EX )
    math(EXPR PACKAGE_CLASSIFICATION_IDX_PLUS_1
      "${PACKAGE_CLASSIFICATION_IDX}+1")
    list(REMOVE_AT ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
      ${PACKAGE_CLASSIFICATION_IDX_PLUS_1} )
  endif()

endmacro()


# @MACRO: tribits_disable_package_on_platforms()
#
# Disable a package automatically for a list of platforms.
#
# Usage::
#
#   tribits_disable_package_on_platforms( <packageName>
#     <hosttype0> <hosttype1> ...)
#
# If any of the host-type arguments ``<hosttypei>`` matches the
# ``${PROJECT_NAME}_HOSTTYPE`` variable for the current platform, then package
# ``<packageName>`` test group classification is changed to ``EX``.  Changing
# the package test group classification to ``EX`` results in the package being
# disabled by default (see `EX packages disabled by default`_).  However,
# an explicit enable can still enable the package.
#
macro( tribits_disable_package_on_platforms  PACKAGE_NAME )
  #message("TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS: ${PACKAGE_NAME}")
  #print_var(${PROJECT_NAME}_HOSTTYPE)
  foreach(HOSTTYPE ${ARGN})
    #print_var(HOSTTYPE)
    if (${PROJECT_NAME}_HOSTTYPE STREQUAL ${HOSTTYPE})
      #message("${${PROJECT_NAME}_HOSTTYPE} == ${HOSTTYPE}")
      tribits_set_package_to_ex(${PACKAGE_NAME})
      #print_var(${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)
      if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
        message(
          "\n***"
          "\n*** NOTE: User has set ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON but the"
          "\n*** package ${PACKAGE_NAME} is not supported on this platform type '${HOSTTYPE}'!"
          "\n***\n"
          )
      endif()
    endif()
  endforeach()
endmacro()


macro( package_disable_on_platforms  PACKAGE_NAME_IN_ )
  message(WARNING "package_disable_on_platforms() is deprecated!"
    "  Use tribits_disable_package_on_platforms() instead!")
  tribits_disable_package_on_platforms(${PACKAGE_NAME_IN_} ${ARGN})
endmacro()



function(tribits_update_ps_pt_ss_st  THING_TYPE  THING_NAME  TESTGROUP_VAR)

  #message("TRIBITS_UPDATE_PS_PT_SS_ST:  ${THING_TYPE}  ${THING_NAME}  ${TESTGROUP_VAR}")

  set(TESTGROUP_IN ${${TESTGROUP_VAR}})
  #print_var(TESTGROUP_IN)

  if (TESTGROUP_IN STREQUAL PS)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "WARNING: ${THING_TYPE} ${THING_NAME} TESTGROUP 'PS' is deprecated."
        "  Use 'PT' instead!")
    endif()
    set(TESTGROUP_OUT PT)
  elseif (TESTGROUP_IN STREQUAL SS)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "WARNING: ${THING_TYPE} ${THING_NAME} TESTGROUP 'SS' is deprecated."
        "  Use 'ST' instead!")
    endif()
    set(TESTGROUP_OUT ST)
  elseif (TESTGROUP_IN STREQUAL TS)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "WARNING: ${THING_TYPE} ${THING_NAME} TESTGROUP 'TS' is deprecated."
        "  Use 'TT' instead!")
    endif()
    set(TESTGROUP_OUT TT)
  else()
    set(TESTGROUP_OUT ${TESTGROUP_IN})
  endif()
  #print_var(TESTGROUP_OUT)

  set(${TESTGROUP_VAR} ${TESTGROUP_OUT} PARENT_SCOPE)

endfunction()
