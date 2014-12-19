# This part needs to be in the top level project CMakeLists.txt.

cmake_minimum_required(VERSION 2.7.20090924)

if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
  message(FATAL_ERROR "Error, CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
    " TrilinosDriver does not support in source builds!" )
endif()

cmake_minimum_required(VERSION 2.7.20090924)

# We have to declare the project name in the top level CMakeLists.txt file.
project(TrilinosDriver NONE)

# The rest of the body is provided by TriBITS.
set(TRIBITS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../tribits)
set(PROJECT_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
include(${TRIBITS_ROOT}/dashboard_driver/TribitsDriverCMakeLists.cmake)

# Provide Trilinos-specific hooks for TRIBITS_DRIVER_SETUP,
# TRIBITS_DRIVER_ADD_DASHBOARD, and
# TRIBITS_ADD_REQUIRED_CMAKE_INSTALLS.
function(TRILINOS_DRIVER_SETUP)
  TRIBITS_DRIVER_SETUP(${ARGN})
endfunction()

function(TRILINOS_DRIVER_ADD_DASHBOARD)
  TRIBITS_DRIVER_ADD_DASHBOARD(${ARGN})
endfunction()

function(TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS)
  TRIBITS_ADD_REQUIRED_CMAKE_INSTALLS(${ARGN})
endfunction()

TDD_PROJECT()
