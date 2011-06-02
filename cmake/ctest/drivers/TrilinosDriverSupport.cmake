include(ParseVariableArguments)
INCLUDE(SetDefaultAndFromEnv)

SET_DEFAULT_AND_FROM_ENV(TDD_FORCE_CMAKE_INSTALL 1)

#
# Placeholder for any global setup necessary on all machines...
#
# Call this function first at the top of any drivers/MACHINE/CMakeLists.txt
# file.
#
function(TRILINOS_DRIVER_SETUP)
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
# TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS after many calls to
# TRILINOS_DRIVER_ADD_DASHBOARD have collected a list of which flavors
# of CMake are required to run all the dashboards on this machine.
#
function(TRILINOS_DRIVER_ADD_TEST_THAT_INSTALLS_CMAKE cmake_type)

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

  if(NOT TD_BASE_DIR OR NOT TRILINOS_HOME_DIR)
    message("error: TD_BASE_DIR and TRILINOS_HOME_DIR must be defined before calling this function")
    return()
  endif()

  find_program(PYTHON_EXE python)

  SET_DEFAULT_AND_FROM_ENV( TDD_HTTP_PROXY "" )
  IF (TDD_HTTP_PROXY)
    SET(TDD_HTTP_PROXY_ARG "--http-proxy=${TDD_HTTP_PROXY}")
  ELSE()
    SET(TDD_HTTP_PROXY_ARG "")
  ENDIF()

  add_test(uninstall-cmake-${cmake_type} ${CMAKE_COMMAND}
    -E remove_directory "${TD_BASE_DIR}/tools/cmake-${cmake_type}"
    )

  add_test(install-cmake-${cmake_type} ${PYTHON_EXE}
    "${TRILINOS_HOME_DIR}/cmake/python/download-cmake.py"
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
#  testname - name of the test as it appears on the TrilinosDriver dashboard
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
function(TRILINOS_DRIVER_ADD_DASHBOARD testname scriptname)

  # Uncomment this line to see output of below PARSE_ARGUMENTS:
  #
  #set(PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED TRUE)

  PARSE_ARGUMENTS(
    # prefix
    PARSE
    # args
    "CTEST_INSTALLER_TYPE;ENVIRONMENT;TIMEOUT_MINUTES;DEPENDS;REQUIRED_FILES"
    # options
    "RUN_SERIAL"
    # the stuff to parse:
    ${ARGN}
  )

  if(NOT "${PARSE_DEFAULT_ARGS}" STREQUAL "")
    message("warning: unrecognized arguments '${PARSE_DEFAULT_ARGS}' to TRILINOS_DRIVER_ADD_DASHBOARD")
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

  #message("TRILINOS_DRIVER_ADD_DASHBOARD")
  #message("  binary_dir='${CMAKE_CURRENT_SOURCE_DIR}'")
  #message("  source_dir='${CMAKE_CURRENT_BINARY_DIR}'")
  #message("  ctest_type='${ctest_type}'")
  #message("  scriptname='${scriptname}'")
  #message("  TD_BASE_DIR='${TD_BASE_DIR}'")
  #message("  testname='${testname}'")
  #message("  timeout_seconds='${timeout_seconds}'")

  add_test(${testname} ${CMAKE_COMMAND}
    -D binary_dir=${CMAKE_CURRENT_BINARY_DIR}
    -D source_dir=${CMAKE_CURRENT_SOURCE_DIR}
    -D ctest_type=${ctest_type}
    -D scriptname=${scriptname}
    -D TD_BASE_DIR=${TD_BASE_DIR}
    -D testname=${testname}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../LocateCTestAndRunScript.cmake
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
    FAIL_REGULAR_EXPRESSION "error: TRILINOS_CTEST_DRIVER_ERROR_QUEUE")

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
# TRILINOS_DRIVER_ADD_DASHBOARD.
#
function(TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS)

  IF (TDD_FORCE_CMAKE_INSTALL STREQUAL 1) 

    get_property(types GLOBAL PROPERTY TD_CMAKE_INSTALLER_TYPES)
  
    if (types)
      list(REMOVE_DUPLICATES types)
      foreach (type ${types})
        MESSAGE(STATUS "Adding CMake install test ${type}")
        TRILINOS_DRIVER_ADD_TEST_THAT_INSTALLS_CMAKE(${type})
      endforeach()
    else()
      message("warning: no cmake install tests are necessary...")
      message("  Put calls to TRILINOS_DRIVER_ADD_DASHBOARD")
      message("  *before* calls to TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS...")
    endif()

  ELSE()

    MESSAGE(STATUS "Skipping CMake install tests because TDD_FORCE_CMAKE_INSTALL!=1")

  ENDIF()

endfunction()
