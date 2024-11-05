# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

#
# CMakeLists.txt body used by the TriBITS driver job.
#

# The base directory is the parent of the binary directory.
# FIXME This is duplicated information. The base directory is set in
# tdd_driver.py. We need to go through two separate CTest invocations
# to get to this point from there. Since environment variables are the
# only way to pass data to CTest, we need to pass this directory
# through the environment.
get_filename_component(TD_BASE_DIR ${CMAKE_BINARY_DIR} PATH)

# CMake versions greater than 2.8.4 have CMAKE_CURRENT_LIST_DIR. Make
# sure it's defined.
if( NOT DEFINED CMAKE_CURRENT_LIST_DIR )
  get_filename_component( CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH )
endif()

# Locate the TriBITS dependencies.
if (NOT TRIBITS_ROOT)
  get_filename_component(TRIBITS_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
endif()
get_filename_component(TRIBITS_ROOT "${TRIBITS_ROOT}" ABSOLUTE)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_LIST_DIR}
  ${TRIBITS_ROOT}/core/utils
  ${TRIBITS_ROOT}/core/package_arch
  ${TRIBITS_ROOT}/core/config_tests
  )

set(TRIBITS_PYTHON_UTILS_DIR "${TRIBITS_ROOT}/python_utils")

include(SetDefaultAndFromEnv)
include(CTest)
include(TribitsDriverSupport)
include(TribitsConfigureCTestCustom)

tribits_configure_ctest_custom(${PROJECT_HOME_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# Function to make exactly one add_subdirectory call based on the site
# name of the machine we're running on. By default, the subdirectory
# is taken to be the site name. The environment variable
# TDD_DRIVER_SUBDIRECTORY can be used to override the default value.
#
function(tdd_project)

  site_name(site)
  set(subdir "${site}")
  message("site='${site}'")

  # But if that directory does not exist as named, and there's a regex match
  # with the name of a subdirectory, use the exact subdirectory name instead:
  #
  file(GLOB filesAndDirs RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*")

  set(dirs "")
  foreach(dir ${filesAndDirs})
    if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}")
      set(dirs ${dirs} "${dir}")
    endif()
  endforeach()

  foreach(dir ${dirs})
    if("${site}" MATCHES "${dir}")
      set(subdir "${dir}")
      message("site='${site}' MATCHES directory name dir='${dir}'")
    endif()
  endforeach()

  # Allow an environment variable to override the test directory.
  set_default_and_from_env(TDD_DRIVER_SUBDIRECTORY "${subdir}")

  # The one add_subdirectory call:
  #
  message("TDD_DRIVER_SUBDIRECTORY='${TDD_DRIVER_SUBDIRECTORY}'")

  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TDD_DRIVER_SUBDIRECTORY}")
    message(FATAL_ERROR "error: there is no subdirectory of ${CMAKE_CURRENT_SOURCE_DIR} matching '${TDD_DRIVER_SUBDIRECTORY}'")
  endif()

  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TDD_DRIVER_SUBDIRECTORY}/CMakeLists.txt")
    message(FATAL_ERROR "error: there is no CMakeLists.txt file in '${CMAKE_CURRENT_SOURCE_DIR}/${TDD_DRIVER_SUBDIRECTORY}'")
  endif()

  add_subdirectory("${TDD_DRIVER_SUBDIRECTORY}")

endfunction()
