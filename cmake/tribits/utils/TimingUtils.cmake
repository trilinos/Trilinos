# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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
# Return the raw time in seconds since epoch, i.e., since 1970-01-01 00:00:00
# UTC
#
FUNCTION(TIMER_GET_RAW_SECONDS   SECONDS_RAW_OUT)
  EXECUTE_PROCESS(COMMAND date "+%s" OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE SECONDS_RAW)
  SET(${SECONDS_RAW_OUT} ${SECONDS_RAW} PARENT_SCOPE)
ENDFUNCTION()


#
# Return the relative time between start and stop seconds
#
FUNCTION(TIMER_GET_REL_SECONDS  SECONDS_RAW_START
  SECONDS_RAW_END  SECONDS_REL_OUT
  )
  MATH(EXPR SECONDS_REL "${SECONDS_RAW_END} - ${SECONDS_RAW_START}")
  SET(${SECONDS_REL_OUT} ${SECONDS_REL} PARENT_SCOPE)
ENDFUNCTION()


#
# Print the relative time in minutes
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
