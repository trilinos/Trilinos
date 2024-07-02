# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(GetCurrentListDir)
include(CMakeParseArguments)
include(SetDefaultAndFromEnv)

# Set if the inner CMake installs are performed or not
set_default_and_from_env(TDD_FORCE_INNER_CMAKE_INSTALL 1)

# Set if the inner CMake installs are performed or not
set_default_and_from_env(TRIBITS_TDD_USE_SYSTEM_CTEST 0)

# Get this value outside of any functions so it will be the path to
# *this* file and not the path to the file calling any of these
# functions.
get_current_list_dir(THIS_SCRIPT_DIR)

#
# Placeholder for any global setup necessary on all machines...
#
# Call this function first at the top of any drivers/MACHINE/CMakeLists.txt
# file.
#
function(tribits_driver_setup)
  #
  # But right now... no global setup required...
  #
endfunction()


#
# Add a test to download and install a single "flavor" of CMake.
#
# Known values of cmake_type are 'min', 'release', 'rc' and 'dev'.
#
# This function should be called automatically by
# TRIBITS_ADD_REQUIRED_CMAKE_INSTALLS after many calls to
# TRIBITS_DRIVER_ADD_DASHBOARD have collected a list of which flavors
# of CMake are required to run all the dashboards on this machine.
#
function(tribits_driver_add_test_that_installs_cmake cmake_type)

  set(known 0)
  if("${cmake_type}" STREQUAL "min" OR
     "${cmake_type}" STREQUAL "release" OR
     "${cmake_type}" STREQUAL "rc" OR
     "${cmake_type}" STREQUAL "dev")
    set(known 1)
  endif()
  if(NOT known)
    message("error: unknown --installer-type '${cmake_type}' -- no test added")
    return()
  endif()

  if(NOT TD_BASE_DIR)
    message(FATAL_ERROR "TD_BASE_DIR must be defined before calling this function")
  endif()

  if(NOT TRIBITS_PYTHON_UTILS_DIR)
    message(FATAL_ERROR "TRIBITS_PYTHON_UTILS_DIR must be defined before calling this function")
  endif()

  find_program(PYTHON_EXE python)

  set_default_and_from_env( TDD_HTTP_PROXY "" )
  if (TDD_HTTP_PROXY)
    set(TDD_HTTP_PROXY_ARG "--http-proxy=${TDD_HTTP_PROXY}")
  else()
    set(TDD_HTTP_PROXY_ARG "")
  endif()

  add_test(uninstall-cmake-${cmake_type} ${CMAKE_COMMAND}
    -E remove_directory "${TD_BASE_DIR}/tools/cmake-${cmake_type}"
    )

  add_test(install-cmake-${cmake_type} ${PYTHON_EXE}
    "${TRIBITS_PYTHON_UTILS_DIR}/download-cmake.py"
    "--skip-detect"
    "--install-dir=${TD_BASE_DIR}/tools/cmake-${cmake_type}"
    "--installer-type=${cmake_type}"
    "${TDD_HTTP_PROXY_ARG}"
    )

  set_property(TEST install-cmake-${cmake_type}
    PROPERTY DEPENDS "uninstall-cmake-${cmake_type}")

endfunction()


#
# Add a test that runs a dashboard script using a known flavor of ctest.
#
# Required arguments:
#  testname - name of the test as it appears on the driver dashboard
#
#  scriptname - name of the file in CMAKE_CURRENT_SOURCE_DIR to run as
#               ctest -S script
#
# Optional arguments:
#
#  [CTEST_INSTALLER_TYPE min|release|rc|dev]
#  [ENVIRONMENT var1=value1;var2=value2;var3=value3;...]
#  [PROCESSORS 1]
#  [RUN_SERIAL]
#  [TIMEOUT_MINUTES 180]
#
function(tribits_driver_add_dashboard testname scriptname)

  message("TRIBITS_DRIVER_ADD_DASHBOARD:  '${testname}'  '${scriptname}' [${ARGN}]")

  cmake_parse_arguments(
    # prefix
    PARSE
    # options
    "RUN_SERIAL"
    # one_value_keywords
    ""
    # multi_value_keywords
    "CTEST_INSTALLER_TYPE;ENVIRONMENT;TIMEOUT_MINUTES;DEPENDS;REQUIRED_FILES"
    # the stuff to parse:
    ${ARGN}
  )

  if(NOT "${PARSE_DEFAULT_ARGS}" STREQUAL "")
    message("warning: unrecognized arguments '${PARSE_DEFAULT_ARGS}' to TRIBITS_DRIVER_ADD_DASHBOARD")
  endif()

  if(PARSE_CTEST_INSTALLER_TYPE)
    set(ctest_type ${PARSE_CTEST_INSTALLER_TYPE})
  else()
    set(ctest_type "release")
  endif()

  if(PARSE_TIMEOUT_MINUTES)
    math(EXPR timeout_seconds "60 * ${PARSE_TIMEOUT_MINUTES}")
  else()
    set(timeout_seconds 10800) # 3 hours...
  endif()

  add_test(${testname} ${CMAKE_COMMAND}
    -D binary_dir=${CMAKE_CURRENT_BINARY_DIR}
    -D source_dir=${CMAKE_CURRENT_SOURCE_DIR}
    -D ctest_type=${ctest_type}
    -D scriptname=${scriptname}
    -D TD_BASE_DIR=${TD_BASE_DIR}
    -D testname=${testname}
    -P ${THIS_SCRIPT_DIR}/LocateCTestAndRunScript.cmake
    )

  # This test that runs the dashboard depends on the test that installs its
  # driving ctest:
  #
  set_property(TEST ${testname} PROPERTY DEPENDS "install-cmake-${ctest_type}")

  if(PARSE_DEPENDS)
    set_property(TEST ${testname} PROPERTY DEPENDS "${PARSE_DEPENDS}")
  endif()

  if(PARSE_REQUIRED_FILES)
    set_property(TEST ${testname} PROPERTY REQUIRED_FILES "${PARSE_REQUIRED_FILES}")
  endif()

  if(PARSE_ENVIRONMENT)
    set_property(TEST ${testname} PROPERTY ENVIRONMENT
      "${PARSE_ENVIRONMENT};TD_CTEST_VERSION_TYPE=${ctest_type}")
  else()
    set_property(TEST ${testname} PROPERTY ENVIRONMENT
      "TD_CTEST_VERSION_TYPE=${ctest_type}")
  endif()

  set_property(TEST ${testname} PROPERTY
    FAIL_REGULAR_EXPRESSION "error: TRIBITS_CTEST_DRIVER_ERROR_QUEUE")

  if(PARSE_PROCESSORS)
    set_property(TEST ${testname} PROPERTY PROCESSORS "${PARSE_PROCESSORS}")
  endif()

  if(PARSE_RUN_SERIAL)
    set_property(TEST ${testname} PROPERTY RUN_SERIAL ON)
  endif()

  set_property(TEST ${testname} PROPERTY TIMEOUT "${timeout_seconds}")

  # Track the required cmake installer types in a global property:
  #
  set_property(GLOBAL APPEND PROPERTY TD_CMAKE_INSTALLER_TYPES "${ctest_type}")

endfunction()


#
# Call this function last at the bottom of any drivers/MACHINE/CMakeLists.txt
# file. It will add tests as needed to download and install all the right
# flavors of ctest required to drive the other tests added with
# TRIBITS_DRIVER_ADD_DASHBOARD.
#
function(tribits_add_required_cmake_installs)

  if (TRIBITS_TDD_USE_SYSTEM_CTEST STREQUAL "1")

    message(STATUS "Skipping CMake install tests because TRIBITS_TDD_USE_SYSTEM_CTEST==1")

  elseif (TDD_FORCE_INNER_CMAKE_INSTALL STREQUAL "1")

    get_property(types GLOBAL PROPERTY TD_CMAKE_INSTALLER_TYPES)

    if (types)
      list(REMOVE_DUPLICATES types)
      foreach (type ${types})
        message(STATUS "Adding CMake install test ${type}")
        tribits_driver_add_test_that_installs_cmake(${type})
      endforeach()
    else()
      message("warning: no cmake install tests are necessary...")
      message("  Put calls to TRIBITS_DRIVER_ADD_DASHBOARD")
      message("  *before* calls to TRIBITS_ADD_REQUIRED_CMAKE_INSTALLS...")
    endif()

  else()

    message(STATUS "Skipping CMake install tests because TDD_FORCE_INNER_CMAKE_INSTALL!=1")

  endif()

endfunction()
