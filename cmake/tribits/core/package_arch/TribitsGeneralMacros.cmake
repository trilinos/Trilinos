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

include(AppendSet)
include(AssertDefined)
include(MessageWrapper)
include(TribitsParseArgumentsHelpers)
include(TribitsSortListAccordingToMasterList)


# Optionally start CMake code configure timing
#
function(tribits_config_code_start_timer  START_TIMER_SECONDS_VAR_OUT)
  if (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    timer_get_raw_seconds(START_TIMER_SECONDS)
    set(${START_TIMER_SECONDS_VAR_OUT} ${START_TIMER_SECONDS} PARENT_SCOPE)
  endif()
endfunction()


# Optionally stop CMake code configure timing
#
function(tribits_config_code_stop_timer  START_TIMER_SECONDS_VAR_IN
  TIMER_STR
  )
  if (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    timer_get_raw_seconds(TIMER_STOP_SECONDS)
    timer_print_rel_time(${${START_TIMER_SECONDS_VAR_IN}}
      ${TIMER_STOP_SECONDS}
      "${TIMER_STR}")
  endif()
endfunction()


# Optionally start CMake code **package** configure timing
#
function(tribits_package_config_code_start_timer  START_TIMER_SECONDS_VAR_OUT)
  if (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING
      AND
      ( ${PROJECT_NAME}_ENABLE_PACKAGE_CONFIGURE_TIMING
        OR ${TRIBITS_PACKAGE}_PACKAGE_CONFIGURE_TIMING )
    )
    timer_get_raw_seconds(START_TIMER_SECONDS)
    set(${START_TIMER_SECONDS_VAR_OUT} ${START_TIMER_SECONDS} PARENT_SCOPE)
  endif()
endfunction()


# Optionally stop CMake code **package** configure timing
#
function(tribits_package_config_code_stop_timer  START_TIMER_SECONDS_VAR_IN
  TIMER_STR
  )
  if (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING
      AND
      ( ${PROJECT_NAME}_ENABLE_PACKAGE_CONFIGURE_TIMING
        OR ${TRIBITS_PACKAGE}_PACKAGE_CONFIGURE_TIMING )
    )
    timer_get_raw_seconds(TIMER_STOP_SECONDS)
    timer_print_rel_time(${${START_TIMER_SECONDS_VAR_IN}}
      ${TIMER_STOP_SECONDS}
      "${TIMER_STR}")
  endif()
endfunction()


# Set the combined directory name taking into account '.' repos.
#
function(tribits_get_repo_name  REPO_DIR  REPO_NAME_OUT)
  if (REPO_DIR STREQUAL ".")
    set(REPO_NAME ${PROJECT_NAME})
  else()
    set(REPO_NAME ${REPO_DIR})
  endif()
  set(${REPO_NAME_OUT} "${REPO_NAME}" PARENT_SCOPE)
endfunction()


# Get the REPO_NAME and REPO_DIR given the REPO
#
function(tribits_get_repo_name_dir  REPO_IN  REPO_NAME_OUT  REPO_DIR_OUT)
  #message("TRIBITS_GET_REPO_NAME_DIR:  '${REPO_IN}'  '${REPO_NAME_OUT}'  '${REPO_DIR_OUT}'")
  # This list of repositories is the list of directories!
  set(REPO_DIR ${REPO_IN})
  # Get the Repository name
  if (REPO_IN STREQUAL ".")
    # The Project and the Reposiotry are one and the same
    set(REPO_NAME ${PROJECT_NAME})
  else()
    # The Repository name is the same as the repository directory
    set(REPO_NAME ${REPO_IN})
  endif()
  set(${REPO_NAME_OUT} ${REPO_NAME} PARENT_SCOPE)
  set(${REPO_DIR_OUT} ${REPO_DIR} PARENT_SCOPE)
endfunction()


# Set the combined directory name taking into account '.' repos.
#
function(tribits_set_base_repo_dir  BASE_DIR  REPO_DIR  BASE_REPO_DIR_OUT)
  if (REPO_DIR STREQUAL ".")
    set(REPO_DIR_STR "")
  else()
    set(REPO_DIR_STR "/${REPO_DIR}")
  endif()
  set(${BASE_REPO_DIR_OUT} "${BASE_DIR}${REPO_DIR_STR}" PARENT_SCOPE)
endfunction()


# Get the list of explicitly enabled entries
#
# These is the list of entries in ${LISTVAR} for which:
#
#   if (${ENABLED_PREFIX}_ENABLE_{ENTRY})
#
# evaluates to true.
#
function(tribits_get_enabled_list  LISTVAR  ENABLED_PREFIX  
  ENABLED_LIST_OUT_OUT  NUM_ENABLED_OUT_OUT
  )
  set(ENABLED_LIST_OUT)
  foreach(ENTITY ${${LISTVAR}})
    set(ENTITY_NAME ${ENABLED_PREFIX}_ENABLE_${ENTITY})
    assert_defined(${ENTITY_NAME})
    set(INCLUDE_ENTITY FALSE)
    if (${ENTITY_NAME})
      list(APPEND  ENABLED_LIST_OUT  ${ENTITY})
    endif()
  endforeach()
  list(LENGTH  ENABLED_LIST_OUT  NUM_ENABLED_OUT)
  set(${ENABLED_LIST_OUT_OUT} ${ENABLED_LIST_OUT} PARENT_SCOPE)
  if (NUM_ENABLED_OUT_OUT)
    set(${NUM_ENABLED_OUT_OUT} ${NUM_ENABLED_OUT} PARENT_SCOPE)
  endif()
endfunction()


# Get the list non-disabled entries
#
# These is the list of entries in ${LISTVAR} for which:
#
#   if (
#     (${ENABLED_PREFIX}_ENABLE_{ENTRY})
#     OR
#     (${ENABLED_PREFIX}_ENABLE_{ENTRY} STREQUAL "" )
#     )
#
# evaluates to true.
#
function(tribits_get_nondisabled_list  LISTVAR  ENABLED_PREFIX  
  NONDISABLED_LIST_OUT_OUT  NUM_NONDISABLED_OUT_OUT
  )
  set(NONDISABLED_LIST_OUT)
  foreach(ENTITY ${${LISTVAR}})
    set(ENTITY_NAME ${ENABLED_PREFIX}_ENABLE_${ENTITY})
    assert_defined(${ENTITY_NAME})
    set(INCLUDE_ENTITY FALSE)
    if (${ENTITY_NAME} OR ${ENTITY_NAME} STREQUAL "")
      list(APPEND  NONDISABLED_LIST_OUT  ${ENTITY})
    endif()
  endforeach()
  list(LENGTH  NONDISABLED_LIST_OUT  NUM_NONDISABLED_OUT)
  set(${NONDISABLED_LIST_OUT_OUT} ${NONDISABLED_LIST_OUT} PARENT_SCOPE)
  if (NUM_NONDISABLED_OUT_OUT)
    set(${NUM_NONDISABLED_OUT_OUT} ${NUM_NONDISABLED_OUT} PARENT_SCOPE)
  endif()
endfunction()


# Get the list of explicitly disabled entries
#
# These is the list of entries in ${LISTVAR} for which:
#
#   if (
#     (NOT ${ENABLED_PREFIX}_ENABLE_{ENTRY})
#     AND
#     (NOT ${ENABLED_PREFIX}_ENABLE_{ENTRY} STREQUAL "" )
#     )
#
# evaluates to true.
#
function(tribits_get_disabled_list  LISTVAR  ENABLED_PREFIX  
  DISABLED_LIST_OUT_OUT  NUM_DISABLED_OUT_OUT
  )
  set(DISABLED_LIST_OUT)
  foreach(ENTITY ${${LISTVAR}})
    set(ENTITY_NAME ${ENABLED_PREFIX}_ENABLE_${ENTITY})
    assert_defined(${ENTITY_NAME})
    set(INCLUDE_ENTITY FALSE)
    if ( (NOT ${ENTITY_NAME}) AND (NOT ${ENTITY_NAME} STREQUAL "") )
      list(APPEND  DISABLED_LIST_OUT  ${ENTITY})
    endif()
  endforeach()
  list(LENGTH  DISABLED_LIST_OUT  NUM_DISABLED_OUT)
  set(${DISABLED_LIST_OUT_OUT} ${DISABLED_LIST_OUT} PARENT_SCOPE)
  if (NUM_DISABLED_OUT_OUT)
    set(${NUM_DISABLED_OUT_OUT} ${NUM_DISABLED_OUT} PARENT_SCOPE)
  endif()
endfunction()


# Get the list of non-enabled entries
#
# These is the list of entries in ${LISTVAR} for which:
#
#   if (NOT ${ENABLED_PREFIX}_ENABLE_{ENTRY})
#
# evaluates to true.
#
function(tribits_get_nonenabled_list  LISTVAR  ENABLED_PREFIX  
  NONENABLED_LIST_OUT_OUT  NUM_NONENABLED_OUT_OUT
  )
  set(NONENABLED_LIST_OUT)
  foreach(ENTITY ${${LISTVAR}})
    set(ENTITY_NAME ${ENABLED_PREFIX}_ENABLE_${ENTITY})
    assert_defined(${ENTITY_NAME})
    set(INCLUDE_ENTITY FALSE)
    if (NOT ${ENTITY_NAME}) # Note that empty "" is also false!
      list(APPEND  NONENABLED_LIST_OUT  ${ENTITY})
    endif()
  endforeach()
  list(LENGTH  NONENABLED_LIST_OUT  NUM_NONENABLED_OUT)
  set(${NONENABLED_LIST_OUT_OUT} ${NONENABLED_LIST_OUT} PARENT_SCOPE)
  if (NUM_NONENABLED_OUT_OUT)
    set(${NUM_NONENABLED_OUT_OUT} ${NUM_NONENABLED_OUT} PARENT_SCOPE)
  endif()
endfunction()


# Macro that sets up the basic lists of enabled packages and packages.
#
macro(tribits_set_up_enabled_lists_and_pkg_idx)

  # ${PROJECT_NAME}_ENABLED_PACKAGES
  tribits_get_enabled_list(
    ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_ENABLED_INTERNAL_TOPLEVEL_PACKAGES
    ${PROJECT_NAME}_NUM_ENABLED_INTERNAL_TOPLEVEL_PACKAGES)

  # ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES
  tribits_get_enabled_list( ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_NUM_ENABLED_INTERNAL_PACKAGES)

  # ${PROJECT_NAME}_REVERSE_ENABLED_INTERNAL_PACKAGES
  set(${PROJECT_NAME}_REVERSE_ENABLED_INTERNAL_PACKAGES
    "${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES}")
  list(REVERSE ${PROJECT_NAME}_REVERSE_ENABLED_INTERNAL_PACKAGES)

  # ${PACKAGE_NAME}_PKG_IDX
  set(PKG_IDX 0)
  foreach(tribitsPackage ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
    set(${tribitsPackage}_PKG_IDX ${PKG_IDX})
    math(EXPR  PKG_IDX  "${PKG_IDX} + 1")
  endforeach()

  # ${PROJECT_NAME}_ENABLED_TPLS
  tribits_get_enabled_list( ${PROJECT_NAME}_DEFINED_TPLS  TPL
    ${PROJECT_NAME}_ENABLED_TPLS  ${PROJECT_NAME}_NUM_ENABLED_TPLS)

  # ${PROJECT_NAME}_REVERSE_ENABLED_TPLS
  set(${PROJECT_NAME}_REVERSE_ENABLED_TPLS
    "${${PROJECT_NAME}_ENABLED_TPLS}")
  list(REVERSE ${PROJECT_NAME}_REVERSE_ENABLED_TPLS)

endmacro()


# @FUNCTION: tribits_set_st_for_dev_mode()
#
# Function that allows packages to easily make a feature ``ST`` for
# development builds and ``PT`` for release builds by default.
#
# Usage::
#
#   tribits_set_st_for_dev_mode(<outputVar>)
#
# This function is typically called in a package's top-level
# `<packageDir>/CMakeLists.txt`_ file before defining other options for the
# package.  The output variable ``${<outputVar>}`` is set to ``ON`` or ``OFF``
# based on the configure state.  In development mode
# (i.e. ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE==ON``), ``${<outputVar>}``
# will be set to ``ON`` only if ``ST`` code is enabled
# (i.e. ``${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE==ON``), otherwise it is
# set to ``OFF``. In release mode
# (i.e. ``${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE==OFF``), ``${<outputVar>}``
# is always set to ``ON``.  This allows some parts of a TriBITS package to be
# considered ``ST`` for development mode (thereby reducing testing time by not
# enabling the dependent features/tests), while still having important
# functionality available to users by default in a release of the package.
#
function(tribits_set_st_for_dev_mode  OUTPUT_VAR)
  if(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
    set(OUTPUT_VAL ${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE})
  else()
    set(OUTPUT_VAL ON)
  endif()
  set(${OUTPUT_VAR} ${OUTPUT_VAL} PARENT_SCOPE)
endfunction()


# For backward compatibility
macro(tribits_set_ss_for_dev_mode  OUTPUT_VAR)
  message(WARNING
    "WARNING: tribits_set_ss_for_dev_mode() is deprecated,"
    " use tribits_set_st_for_dev_mode() instead!")
  tribits_set_st_for_dev_mode(${OUTPUT_VAR})
endmacro()


# Function that extracts all of the required and optional
# items for a given class of package lists
#
function( tribits_gather_enabled_items  PACKAGE_NAME  LISTTYPE_PREFIX
  LISTTYPE_POSTFIX  GATHERED_ITEMS_LIST_OUT
  )

  #message("TRIBITS_GATHER_ENABLED_ITEMS:  '${PACKAGE_NAME}'  '${LISTTYPE_PREFIX}'"
  #  "  '${LISTTYPE_POSTFIX}'  '${GATHERED_ITEMS_LIST_OUT}'")

  set(GATHERED_ITEMS_LIST_TMP
    ${${PACKAGE_NAME}_${LISTTYPE_PREFIX}_REQUIRED_DEP_${LISTTYPE_POSTFIX}}
    )

  #message("TRIBITS_GATHER_ENABLED_ITEMS:"
  #  "  ${PACKAGE_NAME}_${LISTTYPE_PREFIX}_REQUIRED_DEP_${LISTTYPE_POSTFIX} = ${GATHERED_ITEMS_LIST_TMP}")

  foreach(ITEM
    ${${PACKAGE_NAME}_${LISTTYPE_PREFIX}_OPTIONAL_DEP_${LISTTYPE_POSTFIX}}
    )
    assert_defined(${PACKAGE_NAME}_ENABLE_${ITEM})
    if (${PACKAGE_NAME}_ENABLE_${ITEM})
      append_set(GATHERED_ITEMS_LIST_TMP ${ITEM})
    endif()
  endforeach()

  #message("TRIBITS_GATHER_ENABLED_ITEMS:"
  #  "  ${GATHERED_ITEMS_LIST_OUT} = ${GATHERED_ITEMS_LIST_TMP}")

  set(${GATHERED_ITEMS_LIST_OUT} ${GATHERED_ITEMS_LIST_TMP} PARENT_SCOPE)

endfunction()


# @FUNCTION: tribits_trace_file_processing()
#
# Print trace of file processing when
# ``${PROJECT_NAME}_TRACE_FILE_PROCESSING`` is ``TRUE``.
#
# Usage::
#
#   tribits_trace_file_processing( <type> <processingType> <filePath>)
#
# Arguments:
#
# * ``<type>``: Allowed values "PROJECT", "REPOSITORY", "PACKAGE", or "TPL"
#
# * ``<processingType>``: Allowed values "INCLUDE", "ADD_SUBDIR", "READ", or
#   "CONFIGURE"
#
# * ``<filePath>``: Path of the file being processed
#
function(tribits_trace_file_processing  TYPE_IN  PROCESSING_TYPE_IN  FILE_PATH)

  if (${PROJECT_NAME}_TRACE_FILE_PROCESSING)

    if (TYPE_IN STREQUAL "PROJECT")
      set(TYPE_STR "PROJECT   ")
    elseif (TYPE_IN STREQUAL "REPOSITORY")
      set(TYPE_STR "REPOSITORY")
    elseif (TYPE_IN STREQUAL "PACKAGE")
      set(TYPE_STR "PACKAGE   ")
    elseif ("${TYPE_IN}" STREQUAL "TPL")
      set(TYPE_STR "TPL       ")
    else()
      message(FATAL_ERROR
        "Error: TYPE_IN='${TYPE_IN}' is invalid!")
    endif()

    if (PROCESSING_TYPE_IN STREQUAL "INCLUDE")
      set(PROCESSING_TYPE_STR "INCLUDE   ")
    elseif (PROCESSING_TYPE_IN STREQUAL "ADD_SUBDIR")
      set(PROCESSING_TYPE_STR "ADD_SUBDIR")
    elseif (PROCESSING_TYPE_IN STREQUAL "READ")
      set(PROCESSING_TYPE_STR "READ      ")
    elseif (PROCESSING_TYPE_IN STREQUAL "CONFIGURE")
      set(PROCESSING_TYPE_STR "CONFIGURE ")
    else()
      message(FATAL_ERROR
        "Error: PROCESSING_TYPE_IN='${PROCESSING_TYPE_IN}' is invalid!")
    endif()

    message("-- " "File Trace: ${TYPE_STR} ${PROCESSING_TYPE_STR} ${FILE_PATH}")

  endif()

endfunction()
