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


include(TimingUtils)


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
