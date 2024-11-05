# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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
