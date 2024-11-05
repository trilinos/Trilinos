# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include("${CMAKE_CURRENT_LIST_DIR}/../test_support/TribitsAddTest.cmake")

include(TribitsAddExecutable)
include(TribitsDeprecatedHelpers)


#
# Helper Macros
#


macro(tribits_fwd_parse_arg  VAR_TO_SET_OUT  ARGNAME)
  if (PARSE_${ARGNAME})
    set(${VAR_TO_SET_OUT} ${${VAR_TO_SET_OUT}} ${ARGNAME} ${PARSE_${ARGNAME}})
  endif()
endmacro()


macro(tribits_fwd_parse_opt  VAR_TO_SET_OUT  OPTNAME)
  if (PARSE_${OPTNAME})
    set(${VAR_TO_SET_OUT} ${${VAR_TO_SET_OUT}} ${OPTNAME})
  endif()
endmacro()


macro(tribits_add_executable_wrapper)
  if (TRIBITS_ADD_EXECUTABLE_AND_TEST_TEST_MODE)
    set(TRIBITS_ADD_EXECUTABLE_CAPTURE_ARGS ${ARGN} CACHE INTERNAL "")
  else()
    tribits_add_executable(${ARGN})
  endif()
endmacro()


macro(tribits_add_test_wrapper)
  if (TRIBITS_ADD_EXECUTABLE_AND_TEST_TEST_MODE)
    set(TRIBITS_ADD_TEST_CAPTURE_ARGS ${ARGN} CACHE INTERNAL "")
  else()
    tribits_add_test(${ARGN})
  endif()
endmacro()


# @FUNCTION: tribits_add_executable_and_test()
#
# Add an executable and a test (or several tests) all in one shot (just calls
# `tribits_add_executable()`_ followed by `tribits_add_test()`_).
#
# Usage::
#
#   tribits_add_executable_and_test(
#     <exeRootName>  [NOEXEPREFIX]  [NOEXESUFFIX]  [ADD_DIR_TO_NAME]
#     SOURCES <src0> <src1> ...
#     [NAME <testName> | NAME_POSTFIX <testNamePostfix>]
#     [CATEGORIES <category0>  <category1> ...]
#     [HOST <host0> <host1> ...]
#     [XHOST <xhost0> <xhost1> ...]
#     [XHOST_TEST <xhost0> <xhost1> ...]
#     [HOSTTYPE <hosttype0> <hosttype1> ...]
#     [XHOSTTYPE <xhosttype0> <xhosttype1> ...]
#     [XHOSTTYPE_TEST <xhosttype0> <xhosttype1> ...]
#     [EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...]
#     [DISABLED "<messageWhyDisabled>"]
#     [DIRECTORY <dir>]
#     [TESTONLYLIBS <lib0> <lib1> ...]
#     [IMPORTEDLIBS <lib0> <lib1> ...]
#     [COMM [serial] [mpi]]
#     [ARGS "<arg0> <arg1> ..." "<arg2> <arg3> ..." ...]
#     [NUM_MPI_PROCS <numProcs>]
#     [RUN_SERIAL]
#     [LINKER_LANGUAGE (C|CXX|Fortran)]
#     [STANDARD_PASS_OUTPUT
#       | PASS_REGULAR_EXPRESSION "<regex0>;<regex1>;..."]
#     [FAIL_REGULAR_EXPRESSION "<regex0>;<regex1>;..."]
#     [WILL_FAIL]
#     [ENVIRONMENT <var0>=<value0> <var1>=<value1> ...]
#     [INSTALLABLE]
#     [TIMEOUT <maxSeconds>]
#     [LIST_SEPARATOR <sep>]
#     [ADDED_EXE_TARGET_NAME_OUT <exeTargetName>]
#     [ADDED_TESTS_NAMES_OUT <testsNames>]
#     )
#
# This function takes a fairly common set of arguments to
# `tribits_add_executable()`_ and `tribits_add_test()`_ but not the full set
# passed to ``tribits_add_test()``.  See the documentation for
# `tribits_add_executable()`_ and `tribits_add_test()`_ to see which arguments
# are accepted by which functions.
#
# Arguments that are specific to this function and not directly passed on to
# ``tribits_add_executable()`` or ``tribits_add_test()`` include:
#
#   ``XHOST_TEST <xhost0> <xhost1> ...``
#
#     When specified, this disables just running the tests for the named hosts
#     ``<xhost0>``, ``<xhost0>`` etc. but still builds the executable for the
#     test.  These are just passed in through the ``XHOST`` argument to
#     ``tribits_add_test()``.
#
#   ``XHOSTTYPE_TEST <xhosttype0> <hosttype1> ...``
#
#     When specified, this disables just running the tests for the named host
#     types ``<hosttype0>``, ``<hosttype0>``, ..., but still builds the
#     executable for the test.  These are just passed in through the
#     ``XHOSTTYPE`` argument to ``tribits_add_test()``.
#
# This is the function to use for simple test executables that you want to run
# that either takes no arguments or just a simple set of arguments passed in
# through ``ARGS``.  For more flexibility, just use
# ``tribits_add_executable()`` followed by ``tribits_add_test()``.
#
# Finally, the tests are only added if tests are enabled for the package
# (i.e. `${PACKAGE_NAME}_ENABLE_TESTS`_ ``= ON``) and other criteria are met.
# But the test executable will always be added if this function is called,
# regardless of the value of ``${PACKAGE_NAME}_ENABLE_TESTS``.  To avoid
# adding the test (or example) executable when
# ``${PACKAGE_NAME}_ENABLE_TESTS=OFF``, put this command in a subdir under
# ``test/`` or ``example/`` and that subdir with
# `tribits_add_test_directories()`_ or `tribits_add_example_directories()`_,
# respectively.
#
function(tribits_add_executable_and_test EXE_NAME)

  #
  # A) Parse the input arguments
  #

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     "RUN_SERIAL;STANDARD_PASS_OUTPUT;WILL_FAIL;ADD_DIR_TO_NAME;INSTALLABLE;NOEXEPREFIX;NOEXESUFFIX"
     #one_value_keywords
     "DISABLED"
     #mulit_value_keywords
     "SOURCES;DEPLIBS;TESTONLYLIBS;IMPORTEDLIBS;NAME;NAME_POSTFIX;NUM_MPI_PROCS;DIRECTORY;KEYWORDS;COMM;ARGS;PASS_REGULAR_EXPRESSION;FAIL_REGULAR_EXPRESSION;ENVIRONMENT;TIMEOUT;LIST_SEPARATOR;CATEGORIES;HOST;XHOST;XHOST_TEST;HOSTTYPE;XHOSTTYPE;EXCLUDE_IF_NOT_TRUE;XHOSTTYPE_TEST;LINKER_LANGUAGE;TARGET_DEFINES;DEFINES;ADDED_EXE_TARGET_NAME_OUT;ADDED_TESTS_NAMES_OUT"
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("")
    message("TRIBITS_ADD_EXECUTABLE_AND_TEST: ${EXE_NAME} ${ARGN}")
  endif()

  if(PARSE_ADDED_EXE_TARGET_NAME_OUT)
    set(${PARSE_ADDED_EXE_TARGET_NAME_OUT} "" PARENT_SCOPE)
  endif()
  if(PARSE_ADDED_TESTS_NAMES_OUT)
    set(${PARSE_ADDED_TESTS_NAMES_OUT} "" PARENT_SCOPE)
  endif()

  #
  # B) Arguments common to both
  #

  set(COMMON_CALL_ARGS "")
  tribits_fwd_parse_arg(COMMON_CALL_ARGS COMM)
  tribits_fwd_parse_arg(COMMON_CALL_ARGS CATEGORIES)
  tribits_fwd_parse_arg(COMMON_CALL_ARGS HOST)
  tribits_fwd_parse_arg(COMMON_CALL_ARGS XHOST)
  tribits_fwd_parse_arg(COMMON_CALL_ARGS HOSTTYPE)
  tribits_fwd_parse_arg(COMMON_CALL_ARGS XHOSTTYPE)
  tribits_fwd_parse_arg(COMMON_CALL_ARGS EXCLUDE_IF_NOT_TRUE)
  tribits_fwd_parse_opt(COMMON_CALL_ARGS NOEXEPREFIX)
  tribits_fwd_parse_opt(COMMON_CALL_ARGS NOEXESUFFIX)

  #
  # C) tribits_add_executable(...)
  #

  if (PARSE_DEPLIBS)
    tribits_deprecated("DEPLIBS argument of tribits_add_executable_and_test() is deprecated.")
  endif()

  set(CALL_ARGS "")
  tribits_fwd_parse_arg(CALL_ARGS SOURCES)
  tribits_fwd_parse_arg(CALL_ARGS DEPLIBS)  # Deprecated
  tribits_fwd_parse_arg(CALL_ARGS TESTONLYLIBS)
  tribits_fwd_parse_arg(CALL_ARGS IMPORTEDLIBS)
  tribits_fwd_parse_arg(CALL_ARGS DIRECTORY)
  tribits_fwd_parse_opt(CALL_ARGS ADD_DIR_TO_NAME)
  tribits_fwd_parse_arg(CALL_ARGS LINKER_LANGUAGE)
  tribits_fwd_parse_arg(CALL_ARGS TARGET_DEFINES)
  tribits_fwd_parse_arg(CALL_ARGS DEFINES)
  tribits_fwd_parse_opt(CALL_ARGS INSTALLABLE)
  if (PARSE_ADDED_EXE_TARGET_NAME_OUT)
    list(APPEND  CALL_ARGS  ADDED_EXE_TARGET_NAME_OUT  ADDED_EXE_TARGET_NAME)
  endif()

  tribits_add_executable_wrapper(${EXE_NAME} ${COMMON_CALL_ARGS} ${CALL_ARGS})

  if(PARSE_ADDED_EXE_TARGET_NAME_OUT)
    set( ${PARSE_ADDED_EXE_TARGET_NAME_OUT}  "${ADDED_EXE_TARGET_NAME}"
      PARENT_SCOPE )
  endif()

  #
  # D) tribits_add_test(...)
  #

  set(CALL_ARGS "")
  tribits_fwd_parse_arg(CALL_ARGS NAME)
  tribits_fwd_parse_arg(CALL_ARGS NAME_POSTFIX)
  tribits_fwd_parse_arg(CALL_ARGS DIRECTORY)
  tribits_fwd_parse_arg(CALL_ARGS KEYWORDS)
  tribits_fwd_parse_arg(CALL_ARGS NUM_MPI_PROCS)
  tribits_fwd_parse_arg(CALL_ARGS ARGS)
  tribits_fwd_parse_arg(CALL_ARGS PASS_REGULAR_EXPRESSION)
  tribits_fwd_parse_arg(CALL_ARGS FAIL_REGULAR_EXPRESSION)
  tribits_fwd_parse_arg(CALL_ARGS ENVIRONMENT)
  tribits_fwd_parse_arg(CALL_ARGS DISABLED)
  tribits_fwd_parse_opt(CALL_ARGS RUN_SERIAL)
  tribits_fwd_parse_opt(CALL_ARGS STANDARD_PASS_OUTPUT)
  tribits_fwd_parse_opt(CALL_ARGS WILL_FAIL)
  tribits_fwd_parse_arg(CALL_ARGS TIMEOUT)
  tribits_fwd_parse_arg(CALL_ARGS LIST_SEPARATOR)
  tribits_fwd_parse_opt(CALL_ARGS ADD_DIR_TO_NAME)
  tribits_fwd_parse_opt(CALL_ARGS ADDED_TESTS_NAMES_OUT)
  if (PARSE_XHOST_TEST)
    set(CALL_ARGS ${CALL_ARGS} XHOST ${PARSE_XHOST_TEST})
  endif()
  if (PARSE_XHOSTTYPE_TEST)
    set(CALL_ARGS ${CALL_ARGS} XHOSTTYPE ${PARSE_XHOSTTYPE_TEST})
  endif()
  if (PARSE_ADDED_TESTS_NAMES_OUT)
    list(APPEND  CALL_ARGS  ADDED_TESTS_NAMES_OUT  ADDED_TESTS_NAMES)
  endif()

  tribits_add_test_wrapper(${EXE_NAME} ${COMMON_CALL_ARGS} ${CALL_ARGS})

  if(PARSE_ADDED_TESTS_NAMES_OUT)
    set(${PARSE_ADDED_TESTS_NAMES_OUT} "${ADDED_TESTS_NAMES}" PARENT_SCOPE)
  endif()

endfunction()

#  LocalWords:  executables
