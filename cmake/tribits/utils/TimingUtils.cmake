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

#
# @FUNCTION: TIMER_GET_RAW_SECONDS()
# 
# Return the raw time in seconds since epoch, i.e., since 1970-01-01 00:00:00
# UTC.
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
# NOTE: This function runs an external process to run the ``date`` command.
# Therefore, it only works on Unix/Linux and other systems that have a
# standard ``date`` command.  Since this runs an external process, this
# function should only be used to time very course-grained operations
# (i.e. that take longer than a second).
#
FUNCTION(TIMER_GET_RAW_SECONDS   SECONDS_RAW_OUT)
  EXECUTE_PROCESS(COMMAND date "+%s" OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE SECONDS_RAW)
  SET(${SECONDS_RAW_OUT} ${SECONDS_RAW} PARENT_SCOPE)
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
FUNCTION(TIMER_GET_REL_SECONDS  SECONDS_RAW_START
  SECONDS_RAW_END  SECONDS_REL_OUT
  )
  MATH(EXPR SECONDS_REL "${SECONDS_RAW_END} - ${SECONDS_RAW_START}")
  SET(${SECONDS_REL_OUT} ${SECONDS_REL} PARENT_SCOPE)
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
# ``<min>m<sec>s`` format.  This can only resolve times a second or greater
# apart.  If the start and end times are less than a second then ``0m0s`` will
# be printed.
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
#   REAL_EXPENSIVE() time: 0m5s
#
# Again, don't try to time something that takes less than 1 second as it will
# be recorded as ``0m0s``.
#   
FUNCTION(TIMER_PRINT_REL_TIME  SECONDS_RAW_START   SECONDS_RAW_END
  MESSAGE_STR
  )
  TIMER_GET_REL_SECONDS(${SECONDS_RAW_START}  ${SECONDS_RAW_END}
     SECONDS_REL)
  # CMake does not support floating point so I need to do this manually
  MATH(EXPR  MINUTES_REL  "${SECONDS_REL}/60")
  MATH(EXPR  SECONDS_REL_REMAINING  "${SECONDS_REL} - 60*${MINUTES_REL}")
  MESSAGE("${MESSAGE_STR}: ${MINUTES_REL}m${SECONDS_REL_REMAINING}s")
ENDFUNCTION()
