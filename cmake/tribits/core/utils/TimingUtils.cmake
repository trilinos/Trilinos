# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

#
# Timing utilities
#
# Warning: Depends on calling 'date' program so will not be portable to all
# platforms so call with care.
#

include("${CMAKE_CURRENT_LIST_DIR}/Split.cmake")


# @FUNCTION: timer_get_raw_seconds()
#
# Return the raw time in seconds (nano-second accuracy) since epoch, i.e.,
# since 1970-01-01 00:00:00 UTC.
#
# Usage::
#
#   timer_get_raw_seconds(<rawSecondsVar>)
#
# This function is used along with `timer_get_rel_seconds()`_, and
# `timer_print_rel_time()`_ to time big chunks of CMake code for timing and
# profiling purposes.  See `timer_print_rel_time()`_ for more details and an
# example.
#
# NOTE: This function runs an external process with ``execute_process()`` to
# run the ``date`` command.  Therefore, it only works on Unix/Linux and other
# systems that have a standard ``date`` command.  Since this uses
# ``execute_process()``, this function should only be used to time very
# course-grained operations (i.e. that take longer than a second).  If the
# `date` command does not exist, then ${<rawSecondsVar>} will be empty on
# output!
#
function(timer_get_raw_seconds   SECONDS_DEC_OUT)
  execute_process(COMMAND date "+%s.%N" OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE SECONDS_DEC)
  set(${SECONDS_DEC_OUT} ${SECONDS_DEC} PARENT_SCOPE)
endfunction()


# @FUNCTION: timer_get_rel_seconds()
#
# Return the relative time between start and stop seconds.
#
# Usage::
#
#   timer_get_rel_seconds(<startSeconds> <endSeconds> <relSecondsOutVar>)
#
# This simple function computes the relative number of seconds between
# ``<startSeconds>`` and ``<endSeconds>`` (returned from
# `timer_get_raw_seconds()`_) and sets the result in the local variable
# ``<relSecondsOutVar>``.
#
function(timer_get_rel_seconds  ABS_SECONDS_DEC_START
  ABS_SECONDS_DEC_END  SECONDS_REL_OUT
  )
  set(DECIMAL_PLACES 3)
  timer_get_truncated_combined_int_from_decimal_num(
    ${ABS_SECONDS_DEC_START}  ${DECIMAL_PLACES}  SECONDS_START_INT)
  #print_var(SECONDS_START_INT)
  timer_get_truncated_combined_int_from_decimal_num(
    ${ABS_SECONDS_DEC_END}  ${DECIMAL_PLACES}  SECONDS_END_INT)
  #print_var(SECONDS_END_INT)

  math(EXPR SECONDS_REL_INT "${SECONDS_END_INT} - ${SECONDS_START_INT}")

  # ToDo: Check that SECONDS_REL_INT is positive and if not, then raise an
  # error and exit!

  timer_get_decimal_num_frum_truncated_combined_int(
    ${SECONDS_REL_INT}    ${DECIMAL_PLACES}  SECONDS_REL )

  set(${SECONDS_REL_OUT} "${SECONDS_REL}" PARENT_SCOPE)

endfunction()


# @FUNCTION: timer_print_rel_time()
#
# Print the relative time between start and stop timers in ``<min>m<sec>s``
# format.
#
# Usage::
#
#   timer_print_rel_time(<startSeconds> <endSeconds> "<messageStr>")
#
# Differences the raw times ``<startSeconds>`` and ``<endSeconds>``
# (i.e. gotten from `timer_get_raw_seconds()`_) and prints the time in
# ``<min>m<sec>s`` format.
#
# This is meant to be used with `timer_get_raw_seconds()`_ to time expensive
# blocks of CMake code like::
#
#   timer_get_raw_seconds(REAL_EXPENSIVE_START)
#
#   real_expensive(...)
#
#   timer_get_raw_seconds(REAL_EXPENSIVE_END)
#
#   timer_print_rel_time(${REAL_EXPENSIVE_START} ${REAL_EXPENSIVE_END}
#      "real_expensive() time")
#
# This will print something like::
#
#   real_expensive() time: 0m5.235s
#
function(timer_print_rel_time  ABS_SECONDS_DEC_START   ABS_SECONDS_DEC_END  MESSAGE_STR )

  timer_get_rel_seconds(${ABS_SECONDS_DEC_START}  ${ABS_SECONDS_DEC_END}  SECONDS_REL)
  #print_var(SECONDS_REL)

  split(${SECONDS_REL}  "[.]"   SECONDS_REL_ARRAY)
  list(GET  SECONDS_REL_ARRAY  0  SECONDS_REL_S)
  list(GET  SECONDS_REL_ARRAY  1  SECONDS_REL_NS)

  # CMake does not support floating point so I need to do this manually
  math(EXPR  MINUTES_REL  "${SECONDS_REL_S}/60")
  math(EXPR  SECONDS_REL_REMAINING  "${SECONDS_REL_S} - 60*${MINUTES_REL}")
  message("${MESSAGE_STR}: ${MINUTES_REL}m${SECONDS_REL_REMAINING}.${SECONDS_REL_NS}s")
endfunction()
# ToDo: Unit test the above function by calling message_wrapper().


#
# Helper functions
#


#
# Create a combined number of seconds as a combined integer given input
# decimal form <sec>.<nanosec>.  This single int form can be subtracted by
# CMake (the decimal form can't).
#
function(timer_get_truncated_combined_int_from_decimal_num
  ABS_SECONDS_DEC_IN  DECIMAL_PLACES  SECONDS_NS_INT_OUT
  )
  #print_var(ABS_SECONDS_DEC_IN)

  split(${ABS_SECONDS_DEC_IN} "[.]"  ABS_SECONDS_DEC_IN_ARRAY)
  list(GET  ABS_SECONDS_DEC_IN_ARRAY  0  ABS_SECONDS_DEC_IN_S)
  list(GET  ABS_SECONDS_DEC_IN_ARRAY  1  ABS_SECONDS_DEC_IN_NS)
  #print_var(ABS_SECONDS_DEC_IN_S)
  #print_var(ABS_SECONDS_DEC_IN_NS)

  string(SUBSTRING ${ABS_SECONDS_DEC_IN_S} ${DECIMAL_PLACES} -1  SECONDS_S_TRUNCATED)
  string(SUBSTRING ${ABS_SECONDS_DEC_IN_NS} 0 ${DECIMAL_PLACES}  SECONDS_NS_TRUNCATED)
  #print_var(SECONDS_S_TRUNCATED)
  #print_var(SECONDS_NS_TRUNCATED)

  set(${SECONDS_NS_INT_OUT} "${SECONDS_S_TRUNCATED}${SECONDS_NS_TRUNCATED}"
    PARENT_SCOPE)

endfunction()


#
# Reform the formatted seconds in the form <sec>.<nanasec> from the combined
# int form.
#
function(timer_get_decimal_num_frum_truncated_combined_int
  SECONDS_NS_INT_IN  DECIMAL_PLACES
  SECONDS_DEC_IN_OUT
  )
  #print_var(SECONDS_NS_INT_IN)

  string(LENGTH ${SECONDS_NS_INT_IN}  COMBINED_INT_LEN)
  math(EXPR SEC_DIGITS "${COMBINED_INT_LEN}-${DECIMAL_PLACES}")
  #print_var(COMBINED_INT_LEN)
  #print_var(SEC_DIGITS)

  if (SEC_DIGITS GREATER 0)
    # There is at least one digit to the right of the decimal place
    string(SUBSTRING ${SECONDS_NS_INT_IN} 0 ${SEC_DIGITS}  SECONDS_S_TRUNCATED)
    string(SUBSTRING ${SECONDS_NS_INT_IN} ${SEC_DIGITS} -1  SECONDS_NS_TRUNCATED)
  else()
     # There are no digits to the right of th decimal place.  Therefore, there
     # is 0 sec and we need to pad the decimal places with zeros.
     set(SECONDS_S_TRUNCATED 0)

     math(EXPR NUM_ZEROS "${DECIMAL_PLACES}-${COMBINED_INT_LEN}")
     #print_var(NUM_ZEROS)
     set(SECONDS_NS_TRUNCATED ${SECONDS_NS_INT_IN})
     set(LOOP_IDX 0)
     while (LOOP_IDX LESS ${NUM_ZEROS})
       set(SECONDS_NS_TRUNCATED "0${SECONDS_NS_TRUNCATED}")
       math(EXPR LOOP_IDX "${LOOP_IDX}+1")
     endwhile()
  endif()

  #print_var(SECONDS_S_TRUNCATED)
  #print_var(SECONDS_NS_TRUNCATED)

  set(${SECONDS_DEC_IN_OUT} "${SECONDS_S_TRUNCATED}.${SECONDS_NS_TRUNCATED}"
    PARENT_SCOPE)

endfunction()

# ToDo: Combine this basic code with code for scaling TIMEOUT into a simple
# DecimalMath.cmake module.  We need to be able to add, subtract, and multiply
# to begin with.
