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

#
# Timing utilities
#
# Warning: Depends on calling 'date' program so will not be portable to all
# platforms so call with care.
#

INCLUDE(Split)

#
# @FUNCTION: TIMER_GET_RAW_SECONDS()
#
# Return the raw time in seconds (nano-second accuracy) since epoch, i.e.,
# since 1970-01-01 00:00:00 UTC.
#
# Usage::
#
#   TIMER_GET_RAW_SECONDS(<rawSecondsVar>)
#
# This function is used along with `TIMER_GET_REL_SECONDS()`_, and
# `TIMER_PRINT_REL_TIME()`_ to time big chunks of CMake code for timing and
# profiling purposes.  See `TIMER_PRINT_REL_TIME()`_ for more details and an
# example.
#
# NOTE: This function runs an external process with ``EXECUTE_PROCESS()`` to
# run the ``date`` command.  Therefore, it only works on Unix/Linux and other
# systems that have a standard ``date`` command.  Since this uses
# ``EXECUTE_PROCESS()``, this function should only be used to time very
# course-grained operations (i.e. that take longer than a second).  If the
# `date` command does not exist, then ${<rawSecondsVar>} will be empty on
# output!
#
FUNCTION(TIMER_GET_RAW_SECONDS   SECONDS_DEC_OUT)
  EXECUTE_PROCESS(COMMAND date "+%s.%N" OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE SECONDS_DEC)
  SET(${SECONDS_DEC_OUT} ${SECONDS_DEC} PARENT_SCOPE)
ENDFUNCTION()


#
# @FUNCTION: TIMER_GET_REL_SECONDS()
#
# Return the relative time between start and stop seconds.
#
# Usage::
#
#   TIMER_GET_REL_SECONDS(<startSeconds> <endSeconds> <relSecondsOutVar>)
#
# This simple function computes the relative number of seconds between
# ``<startSeconds>`` and ``<endSeconds>`` (returned from
# `TIMER_GET_RAW_SECONDS()`_) and sets the result in the local variable
# ``<relSecondsOutVar>``.
#
FUNCTION(TIMER_GET_REL_SECONDS  ABS_SECONDS_DEC_START
  ABS_SECONDS_DEC_END  SECONDS_REL_OUT
  )
  SET(DECIMAL_PLACES 3)
  TIMER_GET_TRUNCATED_COMBINED_INT_FROM_DECIMAL_NUM(
    ${ABS_SECONDS_DEC_START}  ${DECIMAL_PLACES}  SECONDS_START_INT)
  #PRINT_VAR(SECONDS_START_INT)
  TIMER_GET_TRUNCATED_COMBINED_INT_FROM_DECIMAL_NUM(
    ${ABS_SECONDS_DEC_END}  ${DECIMAL_PLACES}  SECONDS_END_INT)
  #PRINT_VAR(SECONDS_END_INT)

  MATH(EXPR SECONDS_REL_INT "${SECONDS_END_INT} - ${SECONDS_START_INT}")

  # ToDo: Check that SECONDS_REL_INT is positive and if not, then raise an
  # error and exit!

  TIMER_GET_DECIMAL_NUM_FRUM_TRUNCATED_COMBINED_INT(
    ${SECONDS_REL_INT}    ${DECIMAL_PLACES}  SECONDS_REL )

  SET(${SECONDS_REL_OUT} "${SECONDS_REL}" PARENT_SCOPE)

ENDFUNCTION()


#
# @FUNCTION: TIMER_PRINT_REL_TIME()
#
# Print the relative time between start and stop timers in ``<min>m<sec>s``
# format.
#
# Usage::
#
#   TIMER_PRINT_REL_TIME(<startSeconds> <endSeconds> "<messageStr>")
#
# Differences the raw times ``<startSeconds>`` and ``<endSeconds>``
# (i.e. gotten from `TIMER_GET_RAW_SECONDS()`_) and prints the time in
# ``<min>m<sec>s`` format.
#
# This is meant to be used with `TIMER_GET_RAW_SECONDS()`_ to time expensive
# blocks of CMake code like::
#
#   TIMER_GET_RAW_SECONDS(REAL_EXPENSIVE_START)
#
#   REAL_EXPENSIVE(...)
#
#   TIMER_GET_RAW_SECONDS(REAL_EXPENSIVE_END)
#
#   TIMER_PRINT_REL_TIME(${REAL_EXPENSIVE_START} ${REAL_EXPENSIVE_END}
#      "REAL_EXPENSIVE() time")
#
# This will print something like::
#
#   REAL_EXPENSIVE() time: 0m5.235s
#
FUNCTION(TIMER_PRINT_REL_TIME  ABS_SECONDS_DEC_START   ABS_SECONDS_DEC_END  MESSAGE_STR )

  TIMER_GET_REL_SECONDS(${ABS_SECONDS_DEC_START}  ${ABS_SECONDS_DEC_END}  SECONDS_REL)
  #PRINT_VAR(SECONDS_REL)

  SPLIT(${SECONDS_REL}  "[.]"   SECONDS_REL_ARRAY)
  LIST(GET  SECONDS_REL_ARRAY  0  SECONDS_REL_S)
  LIST(GET  SECONDS_REL_ARRAY  1  SECONDS_REL_NS)

  # CMake does not support floating point so I need to do this manually
  MATH(EXPR  MINUTES_REL  "${SECONDS_REL_S}/60")
  MATH(EXPR  SECONDS_REL_REMAINING  "${SECONDS_REL_S} - 60*${MINUTES_REL}")
  MESSAGE("${MESSAGE_STR}: ${MINUTES_REL}m${SECONDS_REL_REMAINING}.${SECONDS_REL_NS}s")
ENDFUNCTION()
# ToDo: Unit test the above function by calling MESSAGE_WRAPPER().


#
# Helper functions
#


#
# Create a combined number of seconds as a combined integer given input
# decimal form <sec>.<nanosec>.  This single int form can be subtracted by
# CMake (the decimal form can't).
#
FUNCTION(TIMER_GET_TRUNCATED_COMBINED_INT_FROM_DECIMAL_NUM
  ABS_SECONDS_DEC_IN  DECIMAL_PLACES  SECONDS_NS_INT_OUT
  )
  #PRINT_VAR(ABS_SECONDS_DEC_IN)

  SPLIT(${ABS_SECONDS_DEC_IN} "[.]"  ABS_SECONDS_DEC_IN_ARRAY)
  LIST(GET  ABS_SECONDS_DEC_IN_ARRAY  0  ABS_SECONDS_DEC_IN_S)
  LIST(GET  ABS_SECONDS_DEC_IN_ARRAY  1  ABS_SECONDS_DEC_IN_NS)
  #PRINT_VAR(ABS_SECONDS_DEC_IN_S)
  #PRINT_VAR(ABS_SECONDS_DEC_IN_NS)

  STRING(SUBSTRING ${ABS_SECONDS_DEC_IN_S} ${DECIMAL_PLACES} -1  SECONDS_S_TRUNCATED)
  STRING(SUBSTRING ${ABS_SECONDS_DEC_IN_NS} 0 ${DECIMAL_PLACES}  SECONDS_NS_TRUNCATED)
  #PRINT_VAR(SECONDS_S_TRUNCATED)
  #PRINT_VAR(SECONDS_NS_TRUNCATED)

  SET(${SECONDS_NS_INT_OUT} "${SECONDS_S_TRUNCATED}${SECONDS_NS_TRUNCATED}"
    PARENT_SCOPE)

ENDFUNCTION()


#
# Reform the formated seconds in the form <sec>.<nanasec> from the combined
# int form.
#
FUNCTION(TIMER_GET_DECIMAL_NUM_FRUM_TRUNCATED_COMBINED_INT
  SECONDS_NS_INT_IN  DECIMAL_PLACES
  SECONDS_DEC_IN_OUT
  )
  #PRINT_VAR(SECONDS_NS_INT_IN)

  STRING(LENGTH ${SECONDS_NS_INT_IN}  COMBINED_INT_LEN)
  MATH(EXPR SEC_DIGITS "${COMBINED_INT_LEN}-${DECIMAL_PLACES}")
  #PRINT_VAR(COMBINED_INT_LEN)
  #PRINT_VAR(SEC_DIGITS)

  IF (SEC_DIGITS GREATER 0)
    # There is at least one digit to the right of the decimal place
    STRING(SUBSTRING ${SECONDS_NS_INT_IN} 0 ${SEC_DIGITS}  SECONDS_S_TRUNCATED)
    STRING(SUBSTRING ${SECONDS_NS_INT_IN} ${SEC_DIGITS} -1  SECONDS_NS_TRUNCATED)
  ELSE()
     # There are no digits to the right of th decimal place.  Therefore, there
     # is 0 sec and we need to pad the decimal places with zeros.
     SET(SECONDS_S_TRUNCATED 0)

     MATH(EXPR NUM_ZEROS "${DECIMAL_PLACES}-${COMBINED_INT_LEN}")
     #PRINT_VAR(NUM_ZEROS)
     SET(SECONDS_NS_TRUNCATED ${SECONDS_NS_INT_IN})
     SET(LOOP_IDX 0)
     WHILE (LOOP_IDX LESS ${NUM_ZEROS})
       SET(SECONDS_NS_TRUNCATED "0${SECONDS_NS_TRUNCATED}")
       MATH(EXPR LOOP_IDX "${LOOP_IDX}+1")
     ENDWHILE()
  ENDIF()

  #PRINT_VAR(SECONDS_S_TRUNCATED)
  #PRINT_VAR(SECONDS_NS_TRUNCATED)

  SET(${SECONDS_DEC_IN_OUT} "${SECONDS_S_TRUNCATED}.${SECONDS_NS_TRUNCATED}"
    PARENT_SCOPE)

ENDFUNCTION()

# ToDo: Combine this basic code with code for scaling TIMEOUT into a simple
# DecimalMath.cmake module.  We need to be able to add, subtract, and multiply
# to begin with.
