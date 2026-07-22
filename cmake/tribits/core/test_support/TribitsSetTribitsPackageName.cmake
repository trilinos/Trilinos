# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include("${CMAKE_CURRENT_LIST_DIR}/../utils/MessageWrapper.cmake")


# Set the TriBITS package name var if it has not already been set
#
macro(tribits_set_tribits_package_name)
  if ("${PACKAGE_NAME}" STREQUAL "")
    if (NOT "${PROJECT_NAME}" STREQUAL "")
      set(PACKAGE_NAME ${PROJECT_NAME})
    else()
       message_wrapper(FATAL_ERROR "Error! Can't set default PACKAGE_NAME because"
	 " PROJECT_NAME is not set!")
    endif()
  endif()
endmacro()
