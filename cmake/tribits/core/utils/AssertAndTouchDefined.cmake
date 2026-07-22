# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

function(assert_and_touch_defined VARS)
  foreach(VAR ${VARS})
    if(NOT DEFINED ${VAR})
      message(SEND_ERROR "Error, the variable ${VAR} is not defined!")
    else()
      # Read the variable so that it will register as being read!
      set(DUMMY_VAR ${${VAR}})
    endif()
  endforeach()
endfunction()
