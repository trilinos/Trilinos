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


INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)


#
# Helper Macros
#


MACRO(TRIBITS_FWD_PARSE_ARG  VAR_TO_SET_OUT  ARGNAME)
  IF (PARSE_${ARGNAME})
    SET(${VAR_TO_SET_OUT} ${${VAR_TO_SET_OUT}} ${ARGNAME} ${PARSE_${ARGNAME}})
  ENDIF()
ENDMACRO()


MACRO(TRIBITS_FWD_PARSE_OPT  VAR_TO_SET_OUT  OPTNAME)
  IF (PARSE_${OPTNAME})
    SET(${VAR_TO_SET_OUT} ${${VAR_TO_SET_OUT}} ${OPTNAME})
  ENDIF()
ENDMACRO()


FUNCTION(TRIBITS_ADD_EXECUTABLE_WRAPPER)
  IF (TRIBITS_ADD_EXECUTABLE_AND_TEST_TEST_MODE)
    SET(TRIBITS_ADD_EXECUTABLE_CAPTURE_ARGS ${ARGN} CACHE INTERNAL "")
  ELSE()
    TRIBITS_ADD_EXECUTABLE(${ARGN})
  ENDIF()
ENDFUNCTION()


FUNCTION(TRIBITS_ADD_TEST_WRAPPER)
  IF (TRIBITS_ADD_EXECUTABLE_AND_TEST_TEST_MODE)
    SET(TRIBITS_ADD_TEST_CAPTURE_ARGS ${ARGN} CACHE INTERNAL "")
  ELSE()
    TRIBITS_ADD_TEST(${ARGN})
  ENDIF()
ENDFUNCTION()

#
# @FUNCTION: TRIBITS_ADD_EXECUTABLE_AND_TEST()
#
# Add an executable and a test (or several tests) all in one shot (just calls
# `TRIBITS_ADD_EXECUTABLE()`_ followed by `TRIBITS_ADD_TEST()`_).
#
# Usage::
#
#   TRIBITS_ADD_EXECUTABLE_AND_TEST(
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
#     [DIRECTORY <dir>]
#     [TESTONLYLIBS <lib0> <lib1> ...]
#     [IMPORTEDLIBS <lib0> <lib1> ...]
#     [COMM [serial] [mpi]]
#     [ARGS "<arg0> <arg1> ..." "<arg2> <arg3> ..." ...]
#     [NUM_MPI_PROCS <numProcs>]
#     [LINKER_LANGUAGE (C|CXX|Fortran)]
#     [STANDARD_PASS_OUTPUT
#       | PASS_REGULAR_EXPRESSION "<regex0>;<regex1>;..."]
#     [FAIL_REGULAR_EXPRESSION "<regex0>;<regex1>;..."]
#     [WILL_FAIL]
#     [ENVIRONMENT <var0>=<value0> <var1>=<value1> ...]
#     [INSTALLABLE]
#     [TIMEOUT <maxSeconds>]
#     )
#
# This function takes a fairly common set of arguments to
# `TRIBITS_ADD_EXECUTABLE()`_ and `TRIBITS_ADD_TEST()`_ but not the full set
# passed to ``TRIBITS_ADD_TEST()``.  See the documentation for
# `TRIBITS_ADD_EXECUTABLE()`_ and `TRIBITS_ADD_TEST()`_ to see which arguments
# are accepted by which functions.
#
# Arguments that are specific to this function and not directly passed on to
# ``TRIBITS_ADD_EXECUTABLE()`` or ``TRIBITS_ADD_TEST()`` include:
#
#   ``XHOST_TEST <xhost0> <xhost1> ...``
#
#     When specified, this disables just running the tests for the named hosts
#     ``<xhost0>``, ``<xhost0>`` etc. but still builds the executable for the
#     test.  These are just passed in through the ``XHOST`` argument to
#     ``TRIBITS_ADD_TEST()``.
#
#   ``XHOSTTYPE_TEST <xhosttype0> <hosttype1> ...``
#
#     When specified, this disables just running the tests for the named host
#     types ``<hosttype0>``, ``<hosttype0>``, ..., but still builds the
#     executable for the test.  These are just passed in through the
#     ``XHOSTTYPE`` argument to ``TRIBITS_ADD_TEST()``.
#
# This is the function to use for simple test executables that you want to run
# that either takes no arguments or just a simple set of arguments passed in
# through ``ARGS``.  For more flexibility, just use
# ``TRIBITS_ADD_EXECUTABLE()`` followed by ``TRIBITS_ADD_TEST()``.
#
FUNCTION(TRIBITS_ADD_EXECUTABLE_AND_TEST EXE_NAME)

  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "SOURCES;DEPLIBS;TESTONLYLIBS;IMPORTEDLIBS;NAME;NAME_POSTFIX;NUM_MPI_PROCS;DIRECTORY;KEYWORDS;COMM;ARGS;NAME;PASS_REGULAR_EXPRESSION;FAIL_REGULAR_EXPRESSION;ENVIRONMENT;CATEGORIES;HOST;XHOST;XHOST_TEST;HOSTTYPE;XHOSTTYPE;XHOSTTYPE_TEST;LINKER_LANGUAGE;DEFINES;ADDED_EXEC_TARGET_NAME_OUT;ADDED_TESTS_NAMES_OUT"
     #options
     "STANDARD_PASS_OUTPUT;WILL_FAIL;TIMEOUT;ADD_DIR_TO_NAME;INSTALLABLE;NOEXEPREFIX;NOEXESUFFIX"
     ${ARGN}
     )

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRIBITS_ADD_EXECUTABLE_AND_TEST: ${EXE_NAME} ${ARGN}")
  ENDIF()

  IF(PARSE_ADDED_EXEC_TARGET_NAME_OUT)
    SET(
      ${PARSE_ADDED_EXEC_TARGET_NAME_OUT}
      PARENT_SCOPE
      )
  ENDIF()
  IF(PARSE_ADDED_TESTS_NAMES_OUT)
    LIST(APPEND
      ${PARSE_ADDED_TESTS_NAMES_OUT}
      ${TEST_NAME_INSTANCE}
      PARENT_SCOPE
      )
  ENDIF()

  #
  # B) Arguments common to both
  #

  SET(COMMON_CALL_ARGS "")
  TRIBITS_FWD_PARSE_ARG(COMMON_CALL_ARGS COMM)
  TRIBITS_FWD_PARSE_ARG(COMMON_CALL_ARGS CATEGORIES)
  TRIBITS_FWD_PARSE_ARG(COMMON_CALL_ARGS HOST)
  TRIBITS_FWD_PARSE_ARG(COMMON_CALL_ARGS XHOST)
  TRIBITS_FWD_PARSE_ARG(COMMON_CALL_ARGS HOSTTYPE)
  TRIBITS_FWD_PARSE_ARG(COMMON_CALL_ARGS XHOSTTYPE)
  TRIBITS_FWD_PARSE_OPT(COMMON_CALL_ARGS NOEXEPREFIX)
  TRIBITS_FWD_PARSE_OPT(COMMON_CALL_ARGS NOEXESUFFIX)

  #
  # C) TribitsAddExecutable(...)
  #

  SET(CALL_ARGS "")
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS SOURCES)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS DEPLIBS)  # Deprecated
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS TESTONLYLIBS)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS IMPORTEDLIBS)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS DIRECTORY)
  TRIBITS_FWD_PARSE_OPT(CALL_ARGS ADD_DIR_TO_NAME)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS LINKER_LANGUAGE)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS DEFINES)
  TRIBITS_FWD_PARSE_OPT(CALL_ARGS INSTALLABLE)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS ADDED_EXEC_TARGET_NAME_OUT)

  TRIBITS_ADD_EXECUTABLE_WRAPPER(${EXE_NAME} ${COMMON_CALL_ARGS} ${CALL_ARGS})

  #
  # D) TribitsAddTest(...)
  #

  SET(CALL_ARGS "")
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS NAME)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS NAME_POSTFIX)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS DIRECTORY)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS KEYWORDS)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS NUM_MPI_PROCS)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS ARGS)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS PASS_REGULAR_EXPRESSION)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS FAIL_REGULAR_EXPRESSION)
  TRIBITS_FWD_PARSE_ARG(CALL_ARGS ENVIRONMENT)
  TRIBITS_FWD_PARSE_OPT(CALL_ARGS STANDARD_PASS_OUTPUT)
  TRIBITS_FWD_PARSE_OPT(CALL_ARGS WILL_FAIL)
  TRIBITS_FWD_PARSE_OPT(CALL_ARGS TIMEOUT)
  TRIBITS_FWD_PARSE_OPT(CALL_ARGS ADD_DIR_TO_NAME)
  TRIBITS_FWD_PARSE_OPT(CALL_ARGS ADDED_TESTS_NAMES_OUT)
  IF (PARSE_XHOST_TEST)
    SET(CALL_ARGS ${CALL_ARGS} XHOST ${PARSE_XHOST_TEST})
  ENDIF()
  IF (PARSE_XHOSTTYPE_TEST)
    SET(CALL_ARGS ${CALL_ARGS} XHOSTTYPE ${PARSE_XHOSTTYPE_TEST})
  ENDIF()

  TRIBITS_ADD_TEST_WRAPPER(${EXE_NAME} ${COMMON_CALL_ARGS} ${CALL_ARGS})

ENDFUNCTION()
