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


include(TribitsAddExecutableTestHelpers)
include(TribitsGeneralMacros)
include(TribitsTestCategories)

include(CMakeParseArguments)
include(GlobalSet)
include(AppendGlobalSet)
include(AppendStringVarWithSep)
include(PrintVar)
include(AdvancedSet)
include(MessageWrapper)
include(TribitsGetCategoriesString)


# Do initialization for test helpers
#
# This must be run just before the packages define their tests and this macro
# must be run in the base-level project scope.
#
macro(tribits_add_test_helpers_init)
  if (TPL_ENABLE_CUDA)
    set(TRIBITS_TEST_EXTRA_ENVIRONMENT CTEST_KOKKOS_DEVICE_TYPE=gpus)
    set(TRIBITS_RESOURCES_PER_PROCESS gpus:1)
  endif()
endmacro()


# Wrapper function for set_tests_properties() to be used in unit testing.
#
function(tribits_set_tests_properties)
  cmake_parse_arguments(PARSE_ARGV 0 FWD "" "" "")
  if (NOT TRIBITS_ADD_TEST_ADD_TEST_UNITTEST)
    set_tests_properties(${FWD_UNPARSED_ARGUMENTS})
  endif()
  if (TRIBITS_SET_TEST_PROPERTIES_CAPTURE_INPUT)
    append_global_set(TRIBITS_SET_TEST_PROPERTIES_INPUT ${FWD_UNPARSED_ARGUMENTS})
  endif()
endfunction()


# Wrapper function for set_property(TEST ...) to be used in unit testing
#
function(tribits_set_test_property)
  cmake_parse_arguments(PARSE_ARGV 0 FWD "" "" "")
  if (NOT TRIBITS_ADD_TEST_ADD_TEST_UNITTEST)
    set_property(TEST ${FWD_UNPARSED_ARGUMENTS})
  endif()
  if (TRIBITS_SET_TEST_PROPERTIES_CAPTURE_INPUT)
    append_global_set(TRIBITS_SET_TEST_PROPERTIES_INPUT ${FWD_UNPARSED_ARGUMENTS})
  endif()
endfunction()


# Scale a timeout by ${PROJECT_NAME}_SCALE_TEST_TIMEOUT
#
# This function will truncate input TIMEOUT_IN but will allow for a simple
# fractional scale factor.  It will also truncate the scale factor to just one
# decimal place.  CMake math(EXPR ...) does not support floating point so I
# have to use just manual integer arithmetic.
#
function(tribits_scale_timeout  TIMEOUT_IN  SCALED_TIMEOUT_OUT)

  if (${PROJECT_NAME}_SCALE_TEST_TIMEOUT
    AND NOT "${${PROJECT_NAME}_SCALE_TEST_TIMEOUT}" EQUAL "1.0"
    )

    # Strip of any fractional part of the timeout
    split("${TIMEOUT_IN}" "[.]" TIMEOUT_IN_ARRAY)
    #print_var(TIMEOUT_IN_ARRAY)
    list(GET TIMEOUT_IN_ARRAY 0 TIMEOUT_IN_TRUNCATED)
    #print_var(TIMEOUT_IN_TRUNCATED)

    # Split the factional part of the scaling factor into SCALE_TEST_INT and
    # SCALE_TEST_INT_FRAC and do the math ourself with integer floating pint
    split("${${PROJECT_NAME}_SCALE_TEST_TIMEOUT}" "[.]" SCALE_TEST_ARRAY)
    list(LENGTH SCALE_TEST_ARRAY SCALE_TEST_ARRAY_LEN)
    #print_var(SCALE_TEST_ARRAY_LEN)

    list(GET SCALE_TEST_ARRAY 0 SCALE_TEST_INT)
    #print_var(SCALE_TEST_INT)

    math(EXPR TIMEOUT_USED
      "${TIMEOUT_IN_TRUNCATED} * ${SCALE_TEST_INT}")
    #print_var(TIMEOUT_USED)

    if ("${SCALE_TEST_ARRAY_LEN}" GREATER 1)
      # Handle the factional part (only take the first digit)
      list(GET SCALE_TEST_ARRAY 1 SCALE_TEST_INT_FRAC_FULL) # Could be more than 1 digit
      #print_var(SCALE_TEST_INT_FRAC_FULL)
      string(SUBSTRING "${SCALE_TEST_INT_FRAC_FULL}" 0 1 SCALE_TEST_INT_FRAC)
      #print_var(SCALE_TEST_INT_FRAC)
      math(EXPR TIMEOUT_USED
        "${TIMEOUT_USED} + (${SCALE_TEST_INT_FRAC} * ${TIMEOUT_IN_TRUNCATED}) / 10")
      #print_var(TIMEOUT_USED)
    endif()

  else()

    set(TIMEOUT_USED "${TIMEOUT_IN}")

  endif()

  set(${SCALED_TIMEOUT_OUT} ${TIMEOUT_USED} PARENT_SCOPE)

endfunction()


# Function that converts a complete string of command-line arguments
# into a form that add_test(...) can correctly deal with.
#
# The main thing this function does is to replace spaces ' ' with
# array separators ';' since this is how add_test(...) expects to deal
# with command-line arguments, as array arguments.  However, this
# function will not do a replacement of ' ' with ';' if a quote is
# active.  This allows you to pass in quoted arguments and have them
# treated as a single argument.
#
function(tribits_convert_cmnd_arg_string_to_add_test_arg_array CMND_ARG_STRING ARG_ARRAY_VARNAME)

  #message("TRIBITS_CONVERT_CMND_ARG_STRING_TO_ADD_TEST_ARG_ARRAY")
  #print_var(CMND_ARG_STRING)
  #print_var(ARG_ARRAY_VARNAME)

  string(LENGTH ${CMND_ARG_STRING} STR_LEN)
  #print_var(STR_LEN)

  math(EXPR STR_LAST_IDX "${STR_LEN}-1")

  set(NEWSTR)

  set(ACTIVE_QUOTE OFF)

  foreach(IDX RANGE ${STR_LAST_IDX})

    string(SUBSTRING ${CMND_ARG_STRING} ${IDX} 1 STR_CHAR)
    #print_var(STR_CHAR)

    if (STR_CHAR STREQUAL "\"")
      if (NOT ACTIVE_QUOTE)
        set(ACTIVE_QUOTE ON)
      else()
        set(ACTIVE_QUOTE OFF)
      endif()
      #print_var(ACTIVE_QUOTE)
    endif()

    if (NOT STR_CHAR STREQUAL " ")
      set(NEWSTR "${NEWSTR}${STR_CHAR}")
    else()
      if (ACTIVE_QUOTE)
        set(NEWSTR "${NEWSTR}${STR_CHAR}")
      else()
        set(NEWSTR "${NEWSTR};")
      endif()
    endif()

  endforeach()

  #print_var(NEWSTR)

  set(${ARG_ARRAY_VARNAME} ${NEWSTR} PARENT_SCOPE)

endfunction()


# Determine if to add the test based on if testing is enabled for the current
# package.
#
function(tribits_add_test_process_enable_tests  ADD_THE_TEST_OUT)
  if(${PACKAGE_NAME}_ENABLE_TESTS)
   set(ADD_THE_TEST TRUE)
  else()
    message_wrapper(
      "-- ${TEST_NAME}: NOT added test because"
      " ${PACKAGE_NAME}_ENABLE_TESTS='${${PACKAGE_NAME}_ENABLE_TESTS}'.")
   set(ADD_THE_TEST FALSE)
  endif()
  set(${ADD_THE_TEST_OUT} ${ADD_THE_TEST} PARENT_SCOPE)
endfunction()


# Determine if to add the test or not based on [X]HOST and [X]HOSTTYPE arguments
#
# Warning: Arguments for [X]HOST and [X]HOSTTYPE arguments are passed in
# implicitly due to scoping of CMake.
#
function(tribits_add_test_process_host_hosttype  ADD_THE_TEST_OUT)

  if ("${${PROJECT_NAME}_HOSTNAME}" STREQUAL "")
    set(${PROJECT_NAME}_HOSTNAME dummy_host)
  endif()

  set(ADD_THE_TEST TRUE)

  if (ADD_THE_TEST)
    if (NOT PARSE_HOST)
      set (PARSE_HOST ${${PROJECT_NAME}_HOSTNAME})
    endif()
    LIST (FIND PARSE_HOST ${${PROJECT_NAME}_HOSTNAME} INDEX_OF_HOSTNAME_IN_HOST_LIST)
    if (${INDEX_OF_HOSTNAME_IN_HOST_LIST} EQUAL -1)
      set(ADD_THE_TEST FALSE)
      set(HOST_MATCH_MSG 
        "-- ${TEST_NAME}: NOT added test because ${PROJECT_NAME}_HOSTNAME='${${PROJECT_NAME}_HOSTNAME}' does not match list HOST='${PARSE_HOST}'!"
        )
    endif()
  endif()

  if (ADD_THE_TEST)
    if (NOT PARSE_XHOST)
      set (PARSE_XHOST NONE)
    endif()
    LIST (FIND PARSE_XHOST ${${PROJECT_NAME}_HOSTNAME} INDEX_OF_HOSTNAME_IN_XHOST_LIST)
    if (NOT ${INDEX_OF_HOSTNAME_IN_XHOST_LIST} EQUAL -1)
      set(ADD_THE_TEST FALSE)
      set(HOST_MATCH_MSG 
        "-- ${TEST_NAME}: NOT added test because ${PROJECT_NAME}_HOSTNAME='${${PROJECT_NAME}_HOSTNAME}' matches list XHOST='${PARSE_XHOST}'!"
        )
    endif()
  endif()

  if (ADD_THE_TEST)
    if (NOT PARSE_HOSTTYPE)
      set(PARSE_HOSTTYPE ${CMAKE_HOST_SYSTEM_NAME})
    endif()
    LIST (FIND PARSE_HOSTTYPE ${CMAKE_HOST_SYSTEM_NAME} INDEX_OF_HOSTSYSTEMNAME_IN_HOSTTYPE_LIST)
    if (${INDEX_OF_HOSTSYSTEMNAME_IN_HOSTTYPE_LIST} EQUAL -1)
      set(ADD_THE_TEST FALSE)
      set(HOST_MATCH_MSG 
        "-- ${TEST_NAME}: NOT added test because CMAKE_HOST_SYSTEM_NAME='${CMAKE_HOST_SYSTEM_NAME}' does not match list HOSTTYPE='${PARSE_HOSTTYPE}'!"
        )
    endif()
  endif()

  if (ADD_THE_TEST)
    if (NOT PARSE_XHOSTTYPE)
      set(PARSE_XHOSTTYPE NONE)
    endif()
    LIST (FIND PARSE_XHOSTTYPE ${CMAKE_HOST_SYSTEM_NAME} INDEX_OF_HOSTSYSTEMNAME_IN_XHOSTTYPE_LIST)
    if (NOT ${INDEX_OF_HOSTSYSTEMNAME_IN_XHOSTTYPE_LIST} EQUAL -1)
      set(ADD_THE_TEST FALSE)
      set(HOST_MATCH_MSG 
        "-- ${TEST_NAME}: NOT added test because CMAKE_HOST_SYSTEM_NAME='${CMAKE_HOST_SYSTEM_NAME}' matches list XHOSTTYPE='${PARSE_XHOSTTYPE}'!"
        )
    endif()
  endif()

  foreach(VAR_NAME ${PARSE_EXCLUDE_IF_NOT_TRUE})
    if (ADD_THE_TEST AND NOT ${VAR_NAME})
      set(ADD_THE_TEST FALSE)
      set(HOST_MATCH_MSG 
        "-- ${TEST_NAME}: NOT added test because EXCLUDE_IF_NOT_TRUE ${VAR_NAME}='${${VAR_NAME}}'!"
      )
    endif()
  endforeach()

  if (HOST_MATCH_MSG AND ${PROJECT_NAME}_TRACE_ADD_TEST)
    message_wrapper("${HOST_MATCH_MSG}")
  endif()

  set(${ADD_THE_TEST_OUT} ${ADD_THE_TEST} PARENT_SCOPE)

endfunction()


# Determine if to add the test or not based on CATEGORIES arguments
#
# Warning: Argument PARSE_CATEGORIES is passed in implicitly due to scoping of
# CMake.
#
function(tribits_add_test_process_categories  ADD_THE_TEST_OUT)

  tribits_filter_and_assert_categories(PARSE_CATEGORIES)
  set(PARSE_CATEGORIES ${PARSE_CATEGORIES} PARENT_SCOPE)

  set(ADD_THE_TEST FALSE)

  # Set the default test-specific category to basic if it is not set
  if (NOT PARSE_CATEGORIES)
    set (PARSE_CATEGORIES BASIC)
  endif()

  #print_var(${PROJECT_NAME}_TEST_CATEGORIES)
  #print_var(PARSE_CATEGORIES)

  if ("${${PROJECT_NAME}_TEST_CATEGORIES}" STREQUAL "")
    # In case this is not a TriBITS project!
    set(${PROJECT_NAME}_TEST_CATEGORIES BASIC)
  endif()

  # Process the test categories
  assert_defined(${PROJECT_NAME}_TEST_CATEGORIES)
  foreach(CATEGORY_USR_SET ${${PROJECT_NAME}_TEST_CATEGORIES})
    #print_var(CATEGORY_USR_SET)
    #print_var(PARSE_CATEGORIES)
    foreach(CATEGORY ${PARSE_CATEGORIES})
      if (CATEGORY STREQUAL ${CATEGORY_USR_SET})
        # Exact match for the category, add the test
        set(ADD_THE_TEST TRUE)
      elseif(CATEGORY STREQUAL "BASIC")
        if (CATEGORY_USR_SET STREQUAL "CONTINUOUS" OR CATEGORY_USR_SET STREQUAL "NIGHTLY"
          OR CATEGORY_USR_SET STREQUAL "HEAVY"
          )
          set(ADD_THE_TEST TRUE)
        endif()
      elseif(CATEGORY STREQUAL "CONTINUOUS")
        if (CATEGORY_USR_SET STREQUAL "NIGHTLY" OR CATEGORY_USR_SET STREQUAL "HEAVY")
          set(ADD_THE_TEST TRUE)
        endif()
      elseif(CATEGORY STREQUAL "NIGHTLY")
        if (CATEGORY_USR_SET STREQUAL "HEAVY")
          set(ADD_THE_TEST TRUE)
        endif()
      else()
        # No matches for the category, don't add the test
      endif()
    endforeach()
  endforeach()

  #print_var(TEST_NAME)
  #print_var(ADD_THE_TEST)
  #print_var(${PROJECT_NAME}_TRACE_ADD_TEST)

  if (TEST_NAME AND NOT ADD_THE_TEST AND ${PROJECT_NAME}_TRACE_ADD_TEST)
    message_wrapper(
      "-- ${TEST_NAME}: NOT added test because ${PROJECT_NAME}_TEST_CATEGORIES='${${PROJECT_NAME}_TEST_CATEGORIES}' does not match this test's CATEGORIES='${PARSE_CATEGORIES}'!"
        )
  endif()


  set(${ADD_THE_TEST_OUT} ${ADD_THE_TEST} PARENT_SCOPE)
  #print_var(${ADD_THE_TEST_OUT})

endfunction()


# FUNCTION: tribits_add_test_get_exe_binary_name()
#
# Get the full name of a package executable given its root name and other
# arguments.
#
# Usage:
#
#   tribits_add_test_get_exe_binary_name(
#     <execRootName>
#     NOEXEPREFIX_IN
#     NOEXESUFFIX_IN
#     ADD_DIR_TO_NAME
#     EXE_BINARY_NAME_OUT
#     )
#
# By default, the full name of the executable is assumed to be::
#
#   ${PACKAGE_NAME}_<execName>${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}
#
function(tribits_add_test_get_exe_binary_name  EXE_NAME_IN
  NOEXEPREFIX_IN  NOEXESUFFIX_IN ADD_DIR_TO_NAME EXE_BINARY_NAME_OUT
  )
  set(EXE_BINARY_NAME "${EXE_NAME_IN}")
  if(PARSE_ADD_DIR_TO_NAME)
    set(DIRECTORY_NAME "")
    tribits_create_name_from_current_source_directory(DIRECTORY_NAME)
    set(EXE_BINARY_NAME "${DIRECTORY_NAME}_${EXE_BINARY_NAME}")
  endif()
  if (NOEXESUFFIX_IN)
    set(EXECUTABLE_SUFFIX "")
  else()
    set(EXECUTABLE_SUFFIX ${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX})
  endif()
  set(EXE_BINARY_NAME "${EXE_BINARY_NAME}${EXECUTABLE_SUFFIX}")
  if(PACKAGE_NAME AND NOT NOEXEPREFIX_IN)
    set(EXE_BINARY_NAME ${PACKAGE_NAME}_${EXE_BINARY_NAME})
  endif()
  set(${EXE_BINARY_NAME_OUT} ${EXE_BINARY_NAME} PARENT_SCOPE)
endfunction()


# Adjust the directory path to an executable for a test
#
function(tribits_add_test_adjust_directory  EXE_BINARY_NAME  DIRECTORY
  EXECUTABLE_PATH_OUT
  )

   set(EXECUTABLE_PATH "${EXE_BINARY_NAME}")

   if (NOT IS_ABSOLUTE ${EXECUTABLE_PATH})

     if (CMAKE_CONFIGURATION_TYPE)
       set(EXECUTABLE_PATH "${CMAKE_CONFIGURATION_TYPE}/${EXECUTABLE_PATH}")
     endif()

     if (DIRECTORY)
       set(EXECUTABLE_PATH "${DIRECTORY}/${EXECUTABLE_PATH}")
     endif()

     if (NOT IS_ABSOLUTE ${EXECUTABLE_PATH})
       set(EXECUTABLE_PATH "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_PATH}")
     endif()

   endif()

  set(${EXECUTABLE_PATH_OUT} ${EXECUTABLE_PATH} PARENT_SCOPE)

endfunction()


# Get the number of MPI processes to use
#
function(tribits_add_test_get_num_procs_used  NUM_MPI_PROCS_IN
  NUM_MPI_PROCS_VAR_NAME_IN  NUM_MPI_PROCS_USED_OUT
  NUM_MPI_PROCS_USED_NAME_OUT
  )
  if (NOT TPL_ENABLE_MPI)
    set(MPI_EXEC_DEFAULT_NUMPROCS 1)
  elseif (NOT DEFINED MPI_EXEC_DEFAULT_NUMPROCS)
    set(MPI_EXEC_DEFAULT_NUMPROCS 1)
  endif()
  if (NOT TPL_ENABLE_MPI)
    set(MPI_EXEC_MAX_NUMPROCS 1)
  elseif (NOT DEFINED MPI_EXEC_MAX_NUMPROCS)
    set(MPI_EXEC_MAX_NUMPROCS 1)
  endif()
  if (NOT NUM_MPI_PROCS_IN)
    set(NUM_MPI_PROCS_IN ${MPI_EXEC_DEFAULT_NUMPROCS})
    set(NUM_PROCS_VAR_NAME "MPI_EXEC_DEFAULT_NUMPROCS") 
  else()
    set(NUM_PROCS_VAR_NAME "${NUM_MPI_PROCS_VAR_NAME_IN}") 
  endif()
  if (${NUM_MPI_PROCS_IN} MATCHES [0-9]+-[0-9]+)
    string(REGEX REPLACE "([0-9]+)-([0-9]+)" "\\1" MIN_NP ${NUM_MPI_PROCS_IN} )
    string(REGEX REPLACE "([0-9]+)-([0-9]+)" "\\2" MAX_NP ${NUM_MPI_PROCS_IN} )
    if(${MIN_NP} LESS ${MPI_EXEC_MAX_NUMPROCS} AND
      ${MAX_NP} GREATER ${MPI_EXEC_MAX_NUMPROCS}
      )
      set(NUM_PROCS_USED ${MPI_EXEC_MAX_NUMPROCS})
    elseif (${MIN_NP} EQUAL ${MPI_EXEC_MAX_NUMPROCS})
      set(NUM_PROCS_USED ${MIN_NP})
    elseif (${MAX_NP} EQUAL ${MPI_EXEC_MAX_NUMPROCS})
      set(NUM_PROCS_USED ${MAX_NP})
    elseif (${MAX_NP} LESS ${MPI_EXEC_MAX_NUMPROCS})
      set(NUM_PROCS_USED ${MAX_NP})
    else()
      # The number of available processors is outside the given range so the
      # test should not be run.
      set(NUM_PROCS_USED -1)
    endif()
  elseif (${NUM_MPI_PROCS_IN} MATCHES [0-9]+,[0-9]+)
    message(SEND_ERROR "The test ${TEST_NAME} can not be added yet"
      " because it we do not yet support the form of"
      " NUM_MPI_PROCS=${NUM_MPI_PROCS_IN}")
      set(NUM_PROCS_USED -1)
  else()
    if(${NUM_MPI_PROCS_IN} GREATER ${MPI_EXEC_MAX_NUMPROCS})
      set(NUM_PROCS_USED -1)
      if (${PROJECT_NAME}_TRACE_ADD_TEST)
        message_wrapper("-- ${TEST_NAME}: NOT added test because ${NUM_PROCS_VAR_NAME}='${NUM_MPI_PROCS_IN}' > MPI_EXEC_MAX_NUMPROCS='${MPI_EXEC_MAX_NUMPROCS}'!")
      endif()
    else()
      set(NUM_PROCS_USED ${NUM_MPI_PROCS_IN})
    endif()
  endif()
  set(${NUM_MPI_PROCS_USED_OUT}  ${NUM_PROCS_USED}  PARENT_SCOPE)
  set(${NUM_MPI_PROCS_USED_NAME_OUT}  ${NUM_PROCS_VAR_NAME}  PARENT_SCOPE)
endfunction()


# Generate the array of arguments for an MPI run
#
# NOTE: The extra test program arguments are passed through ${ARGN}.
#
function( tribits_add_test_get_test_cmnd_array  CMND_ARRAY_OUT
  EXECUTABLE_PATH  NUM_PROCS_USED
  )
  if (TPL_ENABLE_MPI)
    set(${CMND_ARRAY_OUT}
       "${MPI_EXEC}"
       ${MPI_EXEC_PRE_NUMPROCS_FLAGS}
       ${MPI_EXEC_NUMPROCS_FLAG} ${NUM_PROCS_USED}
       ${MPI_EXEC_POST_NUMPROCS_FLAGS}
       "${EXECUTABLE_PATH}"
       ${ARGN}
       PARENT_SCOPE
       )
  else()
    set(${CMND_ARRAY_OUT}
       "${EXECUTABLE_PATH}"
       ${ARGN}
       PARENT_SCOPE
       )
  endif()
endfunction()


# Get the number of cores used by process
#
function(tribits_add_test_get_num_total_cores_used  TEST_NAME_IN
  NUM_TOTAL_CORES_USED_IN  NUM_TOTAL_CORES_USED_NAME_IN
  NUM_PROCS_USED_IN  NUM_PROCS_USED_NAME_IN
  NUM_TOTAL_CORES_USED_OUT  SKIP_TEST_OUT
  )

  set(SKIP_TEST FALSE)

  if (NUM_TOTAL_CORES_USED_IN)

    set(NUM_TOTAL_CORES_USED  ${NUM_TOTAL_CORES_USED_IN})

    if (NUM_TOTAL_CORES_USED_IN  GREATER  MPI_EXEC_MAX_NUMPROCS)
      if (${PROJECT_NAME}_TRACE_ADD_TEST)
        message_wrapper(
          "-- ${TEST_NAME_IN}: NOT added test because ${NUM_TOTAL_CORES_USED_NAME_IN}='${NUM_TOTAL_CORES_USED_IN}' > MPI_EXEC_MAX_NUMPROCS='${MPI_EXEC_MAX_NUMPROCS}'!"
          )
      endif()
      set(SKIP_TEST TRUE)
    endif()

    if (NUM_PROCS_USED_IN  GREATER  NUM_TOTAL_CORES_USED_IN)
      message_wrapper(
        FATAL_ERROR
        "ERROR: ${TEST_NAME_IN}: ${NUM_PROCS_USED_NAME_IN}='${NUM_PROCS_USED}' > ${NUM_TOTAL_CORES_USED_NAME_IN}='${NUM_TOTAL_CORES_USED_IN}' not allowed!"
        )
      set(SKIP_TEST TRUE) # For unit testing since above will not abort!
    endif()

  else()

    set(NUM_TOTAL_CORES_USED  ${NUM_PROCS_USED_IN})

  endif()

  set(${NUM_TOTAL_CORES_USED_OUT}  ${NUM_TOTAL_CORES_USED}  PARENT_SCOPE)
  set(${SKIP_TEST_OUT}  ${SKIP_TEST}  PARENT_SCOPE)

endfunction()


# Read ${TEST_NAME_IN}_SET_DISABLED_AND_MSG and set output var
# <setDisabledAndMsgOut> as used by the TRIBITS_ADD[_ADVANCED]_test()
# functions.
#
# Usage:
#
#    tribits_set_disabled_and_msg(
#      <testName>              # The full test name passed to add_test()
#      "${PARSE_DISABLED}"     # From the input arg "DISABLED <msg>" to TA[A]t()
#      <setDisabledAndMsgOut>  # Sets var of this name on output
#      )
#
function(tribits_set_disabled_and_msg  TEST_NAME_IN  PARSE_DISABLED
  SET_DISABLED_AND_MSG_OUT
  )

  set(SET_DISABLED_AND_MSG "${PARSE_DISABLED}")
  if (NOT "${${TEST_NAME_IN}_SET_DISABLED_AND_MSG}" STREQUAL "") # Override!
    set(SET_DISABLED_AND_MSG "${${TEST_NAME_IN}_SET_DISABLED_AND_MSG}")
  endif()

  set(${SET_DISABLED_AND_MSG_OUT} "${SET_DISABLED_AND_MSG}" PARENT_SCOPE)

endfunction()


# Read ${TEST_NAME_IN}_SET_RUN_SERIAL and set output var SET_RUN_SERIAL_OUT as
# used by the TRIBITS_ADD[_ADVANCED]_test() functions.
#
# Usage:
#
#    tribits_set_run_serial(
#      <testName>              # The full test name passed to add_test()
#      "${PARSE_RUN_SERIAL}"   # From the input option "RUN_SERIAL" to TA[A]t()
#      <setRunSerialOut>       # Sets var of this name on output
#      )
#
function(tribits_set_run_serial  TEST_NAME_IN  PARSE_RUN_SERIAL
  SET_RUN_SERIAL_OUT
  )

  set(SET_RUN_SERIAL "${PARSE_RUN_SERIAL}")
  if (NOT "${${TEST_NAME_IN}_SET_RUN_SERIAL}" STREQUAL "") # Override!
    set(SET_RUN_SERIAL "${${TEST_NAME_IN}_SET_RUN_SERIAL}")
  endif()

  set(${SET_RUN_SERIAL_OUT} "${SET_RUN_SERIAL}" PARENT_SCOPE)

endfunction()


# Determine if the test should be skipped due to a disable var set
#
# Usage:
#
#   tribits_add_test_query_disable(
#     <disableThisTestVarOut>  # Set var of this name on output
#     )
#
function(tribits_add_test_query_disable  DISABLE_TEST_VAR_OUT)

  #message("tribits_add_test_query_disable(): ${DISABLE_TEST_VAR_OUT} ${ARGN}")
  list(GET ARGN 0 TEST_NAME_IN)
  set(TEST_NAME_DISABLE_VAR  ${TEST_NAME_IN}_DISABLE)
  #print_var(${TEST_NAME_DISABLE_VAR})
  if (${TEST_NAME_DISABLE_VAR})
    if (${PROJECT_NAME}_TRACE_ADD_TEST)
      message_wrapper(
        "-- ${TEST_NAME_IN}: NOT added test because ${TEST_NAME_DISABLE_VAR}='${${TEST_NAME_DISABLE_VAR}}'!")
    endif()
    set(${DISABLE_TEST_VAR_OUT} ON PARENT_SCOPE)
  else()
    set(${DISABLE_TEST_VAR_OUT} OFF PARENT_SCOPE)
  endif()

endfunction()


#
# Determine if to disable test due to <Package>_SKIP_CTEST_ADD_TEST=ON
#
function(tribits_add_test_process_skip_ctest_add_test  ADD_THE_TEST_OUT)
  if(${PACKAGE_NAME}_SKIP_CTEST_ADD_TEST OR ${PARENT_PACKAGE_NAME}_SKIP_CTEST_ADD_TEST)
    if (PARENT_PACKAGE_NAME STREQUAL PACKAGE_NAME)
      set(DISABLE_VAR_MSG
        "${PACKAGE_NAME}_SKIP_CTEST_ADD_TEST='${${PACKAGE_NAME}_SKIP_CTEST_ADD_TEST}'")
    else()
      set(DISABLE_VAR_MSG
        "${PARENT_PACKAGE_NAME}_SKIP_CTEST_ADD_TEST='${${PARENT_PACKAGE_NAME}_SKIP_CTEST_ADD_TEST}'")
    endif()
    message_wrapper(
      "-- ${TEST_NAME}: NOT added test because ${DISABLE_VAR_MSG}!")
   set(ADD_THE_TEST FALSE)
  else()
   set(ADD_THE_TEST TRUE)
  endif()
  set(${ADD_THE_TEST_OUT} ${ADD_THE_TEST} PARENT_SCOPE)
endfunction()


#
# Wrapper for adding a test to facilitate unit testing
#
function(tribits_add_test_add_test TEST_NAME EXE_NAME)

  string(REPLACE "${PARSE_LIST_SEPARATOR}" "\\;" args "${ARGN}")

  if (TRIBITS_ADD_TEST_ADD_TEST_CAPTURE)
    append_global_set(TRIBITS_ADD_TEST_ADD_TEST_INPUT NAME ${TEST_NAME} COMMAND ${args})
  endif()

  if (NOT TRIBITS_ADD_TEST_ADD_TEST_UNITTEST)
    add_test(NAME ${TEST_NAME} COMMAND ${args})
  endif()

  tribits_set_tests_properties(${TEST_NAME} PROPERTIES REQUIRED_FILES ${EXE_NAME})

endfunction()


#
# Set the pass/fail properties of a test that has already been added
#

function(tribits_private_add_test_set_passfail_properties TEST_NAME_IN)

  # PASS_REGULAR_EXPRESSION

  if (PARSE_STANDARD_PASS_OUTPUT)
    tribits_set_tests_properties(${TEST_NAME_IN} PROPERTIES PASS_REGULAR_EXPRESSION
      "End Result: TEST PASSED")
  endif()

  if (PARSE_PASS_REGULAR_EXPRESSION)
    tribits_set_tests_properties(${TEST_NAME_IN} PROPERTIES PASS_REGULAR_EXPRESSION
      "${PARSE_PASS_REGULAR_EXPRESSION}")
  endif()

  if (PARSE_WILL_FAIL)
    tribits_set_tests_properties(${TEST_NAME_IN} PROPERTIES WILL_FAIL ON)
  endif()

  # FAIL_REGULAR_EXPRESSION

  if (PARSE_FAIL_REGULAR_EXPRESSION)
    tribits_set_test_property(${TEST_NAME_IN} APPEND PROPERTY FAIL_REGULAR_EXPRESSION
      "${PARSE_FAIL_REGULAR_EXPRESSION}")
  endif()

  if (${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OR
    ${PROJECT_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE
    )
    tribits_set_test_property(${TEST_NAME_IN} APPEND PROPERTY FAIL_REGULAR_EXPRESSION
      "The following Teuchos::RCPNode objects were created")
    # NOTE: The above string must be kept in sync with the C++ code!
  endif()
  # ToDo: Make a variable ${PROJECT_NAME}_EXTRA_FAIL_REGULAR_EXPRESSION and
  # move the above logic to Trilinos somehow.

endfunction()


#
# Set the timeout for a test already added
#
function(tribits_private_add_test_set_timeout  TEST_NAME_IN   TIMEOUT_USED_OUT)

  if (PARSE_TIMEOUT)
    tribits_scale_timeout("${PARSE_TIMEOUT}" TIMEOUT_USED)
    tribits_set_tests_properties(${TEST_NAME_IN} PROPERTIES TIMEOUT ${TIMEOUT_USED})
  else()
    set(TIMEOUT_USED "")
  endif()

  set(${TIMEOUT_USED_OUT} ${TIMEOUT_USED} PARENT_SCOPE)

endfunction()


#
# Set the environment for a test already added
#
function(tribits_private_add_test_set_environment  TEST_NAME_IN)

  if (PARSE_ENVIRONMENT)
    string(REPLACE "${PARSE_LIST_SEPARATOR}" "\\;" PARSE_ENVIRONMENT
      "${PARSE_ENVIRONMENT}")
    tribits_set_test_property(${TEST_NAME_IN} PROPERTY ENVIRONMENT
      "${PARSE_ENVIRONMENT}")
  endif()

endfunction()


#
# Set the environment for a test already added
#
function(tribits_private_add_test_set_processors  TEST_NAME_IN
  NUM_TOTAL_CORES_USED_IN  PROCESSORS_OUT
  )

  set(PROCESSORS_USED ${NUM_TOTAL_CORES_USED_IN})

  tribits_set_tests_properties(${TEST_NAME_IN} PROPERTIES PROCESSORS
    "${PROCESSORS_USED}")

  set(${PROCESSORS_OUT} ${PROCESSORS_USED} PARENT_SCOPE)

endfunction()


# Print test added message!
#
function(tribits_private_add_test_print_added  TEST_NAME_IN  CATEGORIES_IN
  NUM_MPI_PROCS_IN  PROCESSORS_IN  TIMEOUT_IN  RUN_SERIAL_IN  DISABLED_MSG_IN
  )

  set(ADDED_TEST_PROPS "")
  tribits_get_categories_string("${CATEGORIES_IN}" CATEGORIES_IN_COMMAS)
  append_string_var_with_sep(ADDED_TEST_PROPS
     ", "  "${CATEGORIES_IN_COMMAS}")
  if (NUM_MPI_PROCS_IN AND TPL_ENABLE_MPI)
    append_string_var_with_sep(ADDED_TEST_PROPS
      ", "  "NUM_MPI_PROCS=${NUM_MPI_PROCS_IN}")
  endif()
  if (PROCESSORS_IN)
    append_string_var_with_sep(ADDED_TEST_PROPS
      ", "  "PROCESSORS=${PROCESSORS_IN}")
  endif()
  if (TIMEOUT_IN)
    append_string_var_with_sep(ADDED_TEST_PROPS
      ", "  "TIMEOUT=${TIMEOUT_IN}")
  endif()

  if (RUN_SERIAL_IN)
    append_string_var_with_sep(ADDED_TEST_PROPS
      ", "  "RUN_SERIAL")
  endif()

  if (DISABLED_MSG_IN)
    append_string_var_with_sep(ADDED_TEST_PROPS
      ", "  "DISABLED")
  endif()

  if (ADDED_TEST_PROPS)
   set(ADDED_TEST_PROPS " (${ADDED_TEST_PROPS})") 
  endif()

  if (${PROJECT_NAME}_TRACE_ADD_TEST)
    message_wrapper("-- ${TEST_NAME_IN}: Added test${ADDED_TEST_PROPS}!")
    if (DISABLED_MSG_IN)
     message_wrapper("--  => Reason DISABLED: ${DISABLED_MSG_IN}")
    endif()
  endif()

endfunction()


# Overall add test command
#
# NOTE: Pass the command arguments on the end in ARGSN.
#
function(tribits_add_test_add_test_all  TEST_NAME_IN
  EXECUTABLE_PATH_IN  CATEGORIES_IN  NUM_PROCS_USED_IN  NUM_TOTAL_CORES_USED_IN
  RUN_SERIAL_IN  DISABLED_MSG_IN
  ADDED_TEST_NAME_OUT
  )

  tribits_add_test_get_test_cmnd_array( CMND_ARRAY
    "${EXECUTABLE_PATH_IN}"  "${NUM_PROCS_USED_IN}" ${ARGN} )

  tribits_add_test_query_disable(DISABLE_THIS_TEST  ${TEST_NAME_IN})

  if (NOT  DISABLE_THIS_TEST)

    tribits_add_test_add_test(${TEST_NAME_IN}  ${EXECUTABLE_PATH_IN}  ${CMND_ARRAY})
    set(${ADDED_TEST_NAME_OUT}  ${TEST_NAME_IN}  PARENT_SCOPE)

    tribits_private_add_test_post_process_added_test(${TEST_NAME_IN}
      "${CATEGORIES_IN}" ${NUM_PROCS_USED_IN}  "${NUM_TOTAL_CORES_USED_IN}"
      "${RUN_SERIAL_IN}" "${DISABLED_MSG_IN}")

  else()

    set(${ADDED_TEST_NAME_OUT} "" PARENT_SCOPE)

  endif()

endfunction()


# Set the label and keywords
#
function(tribits_private_add_test_add_label_and_keywords  TEST_NAME_IN)

  tribits_set_test_property(${TEST_NAME_IN} APPEND PROPERTY
    LABELS ${PARENT_PACKAGE_NAME})

  if(PARSE_KEYWORDS)
    tribits_set_test_property(${TEST_NAME_IN} APPEND PROPERTY
      LABELS ${PARSE_KEYWORDS})
  endif()

endfunction()


# Postprocess a test that was added
#
function(tribits_private_add_test_post_process_added_test  TEST_NAME_IN
  CATEGORIES_IN  NUM_PROCS_USED_IN  NUM_TOTAL_CORES_USED_IN
  RUN_SERIAL_IN  DISABLED_MSG_IN
  )

  tribits_private_add_test_set_passfail_properties(${TEST_NAME_IN})
  tribits_private_add_test_set_environment(${TEST_NAME_IN})
  tribits_private_add_test_set_processors(${TEST_NAME_IN}
    "${NUM_TOTAL_CORES_USED_IN}"  PROCESSORS_USED)
  tribits_private_add_test_set_timeout(${TEST_NAME_IN}  TIMEOUT_USED)

  if(RUN_SERIAL_IN)
    tribits_set_tests_properties(${TEST_NAME_IN} PROPERTIES RUN_SERIAL ON)
  endif()

  if(DISABLED_MSG_IN)
    tribits_set_tests_properties(${TEST_NAME_IN} PROPERTIES DISABLED ON)
  endif()

  tribits_private_add_test_add_environment_and_resource(${TEST_NAME_IN}
     ${NUM_PROCS_USED_IN})

  tribits_private_add_test_add_label_and_keywords(${TEST_NAME_IN})

  tribits_private_add_test_print_added(${TEST_NAME_IN}
    "${CATEGORIES_IN}"  "${NUM_PROCS_USED_IN}"  "${PROCESSORS_USED}"
    "${TIMEOUT_USED}" "${RUN_SERIAL_IN}" "${DISABLED_MSG_IN}")

endfunction()


# Add environment and resource properties to a test
#
function(tribits_private_add_test_add_environment_and_resource  TEST_NAME_IN
  NUM_PROCS_USED_IN
  )
  if(TRIBITS_TEST_EXTRA_ENVIRONMENT)
    tribits_set_test_property(${TEST_NAME_IN} APPEND PROPERTY ENVIRONMENT
      "${TRIBITS_TEST_EXTRA_ENVIRONMENT}")
  endif()

  if(TRIBITS_RESOURCES_PER_PROCESS)
    set(NUM_PROCESSES ${NUM_PROCS_USED_IN})
    if(NOT NUM_PROCESSES OR NUM_PROCESSES LESS 1)
      set(NUM_PROCESSES 1)
    endif()
    tribits_set_tests_properties(${TEST_NAME_IN} PROPERTIES RESOURCE_GROUPS
      "${NUM_PROCESSES},${TRIBITS_RESOURCES_PER_PROCESS}")
  endif()
endfunction()
