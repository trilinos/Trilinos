# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Top-level include file that pulls in TriBITS so it can be used by a project.
# A Project's top-level CMakeLists.cmake file just does:
#
#   include(${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake)
#
# and then they can call:
#
#    tribits_project()

if (${PROJECT_NAME}_TRIBITS_DIR)
  set(TRIBITS_BASE_DIR_LOCAL "${${PROJECT_NAME}_TRIBITS_DIR}")
elseif(CMAKE_CURRENT_LIST_DIR)
  set(TRIBITS_BASE_DIR_LOCAL "${CMAKE_CURRENT_LIST_DIR}")
else()
  message(FATAL_ERROR "Please set ${PROJECT_NAME}_TRIBITS_DIR!")
endif()

include("${TRIBITS_BASE_DIR_LOCAL}/core/package_arch/TribitsProject.cmake")
include("${TRIBITS_BASE_DIR_LOCAL}/Version.cmake")
