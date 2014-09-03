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


INCLUDE(TribitsAddAdvancedTestHelpers)

INCLUDE(AppendStringVar)
INCLUDE(PrintVar)


#
# @FUNCTION: TRIBITS_ADD_ADVANCED_TEST()
# 
# Function that creates an advanced test defined by stringing together one or
# more executables and/or commands that is run as a ``cmake -P`` script with
# very flexible pass/fail criteria.
#
# Usage::
#
#   TRIBITS_ADD_ADVANCED_TEST(
#     <testName>
#     TEST_0 (EXEC <execTarget0> | CMND <cmndExec0>) ...
#     [TEST_1 (EXEC <execTarget1> | CMND <cmndExec1>) ...]
#     ...
#     [TEST_N (EXEC <execTargetN> | CMND <cmndExecN>) ...]
#     [OVERALL_WORKING_DIRECTORY (<overallWorkingDir> | TEST_NAME)]
#     [FAIL_FAST]
#     [KEYWORDS <keyword1> <keyword2> ...]
#     [COMM [serial] [mpi]]
#     [OVERALL_NUM_MPI_PROCS <overallNumProcs>]
#     [CATEGORIES <category0> <category1> ...]
#     [HOST <host0> <host1> ...]
#     [XHOST <host0> <host1> ...]
#     [HOSTTYPE <hosttype0> <hosttype1> ...]
#     [XHOSTTYPE <hosttype0> <hosttype1> ...]
#     [FINAL_PASS_REGULAR_EXPRESSION <regex> |
#       FINAL_FAIL_REGULAR_EXPRESSION <regex>]
#     [ENVIRONMENT <var1>=<value1> <var2>=<value2> ...]
#     [TIMEOUT <maxSeconds>]
#     )
#
# This function allows one to add a single CTest test that is actually a
# sequence of one or more separate commands strung together in some way to
# define the final pass/fail. One will want to use this function to add a test
# instead of `TRIBITS_ADD_TEST()`_ when one needs to run more than one
# command, or one needs more sophisticated checking of the test result other
# than just grepping STDOUT (e.g. by running separate post-processing programs
# to examine output files).
#
# Each atomic test case is either a package-built executable or just a basic
# command.  An atomic test command block ``TEST_<idx>`` (i.e. ``TEST_0``,
# ``TEST_1``, ...) takes the form::
#
#   TEST_<idx>
#      (EXEC <exeRootName> [NOEXEPREFIX] [NOEXESUFFIX] [ADD_DIR_TO_NAME]
#             [DIRECTORY <dir>]
#         | CMND <cmndExec>)
#      [ARGS <arg1> <arg2> ... <argn>]
#      [MESSAGE "<message>"]
#      [WORKING_DIRECTORY <workingDir>]
#      [NUM_MPI_PROCS <numProcs>]
#      [OUTPUT_FILE <outputFile>]
#      [NO_ECHO_OUTPUT]]
#      [PASS_ANY
#        | PASS_REGULAR_EXPRESSION "<regex>"
#        | PASS_REGULAR_EXPRESSION_ALL "<regex1>" "<regex2>" ... "<regexn>"
#        | FAIL_REGULAR_EXPRESSION "<regex>"
#        | STANDARD_PASS_OUTPUT
#        ]
#
# By default, each and every atomic test or command needs to pass (as defined below) in
# order for the overall test to pass.
#
# *Sections:*
# 
# * `Overall Arguments (TRIBITS_ADD_ADVANCED_TEST())`_
# * `TEST_<idx> Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Overall Pass/Fail (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Argument Parsing and Ordering (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Implementation Details (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Setting Additional Test Properties (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Running multiple tests at the same time (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Disabling Tests Externally (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Debugging and Examining Test Generation (TRIBITS_ADD_ADVANCED_TEST())`_
#
# .. _Overall Arguments (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Overall Arguments (TRIBITS_ADD_ADVANCED_TEST())**
#
# Below, some of the overall arguments are described.  The rest of the overall
# arguments that control overall pass/fail are described in `Overall Pass/Fail
# (TRIBITS_ADD_ADVANCED_TEST())`_.  (NOTE: All of these arguments must be
# listed outside of the ``TEST_<idx>`` blocks, see `Argument Parsing and
# Ordering (TRIBITS_ADD_ADVANCED_TEST())`_).
#
#   ``<testName>``
#
#     The name of the test (which will have ``${PACKAGE_NAME}_`` prepended to
#     the name) that will be used to name the output CMake script file as well
#     as the CTest test name passed into ``ADD_TEST()``.  This must be the
#     first argument to this function.
#
#   ``OVERALL_WORKING_DIRECTORY <overallWorkingDir>``
#
#     If specified, then the working directory ``<overallWorkingDir>`` will be
#     created and all of the test commands by default will be run from within
#     this directory.  If the value ``<overallWorkingDir>=TEST_NAME`` is
#     given, then the working directory will be given the name
#     ``${PACKAGE_NAME}_<testName>``.  If the directory
#     ``<overallWorkingDir>`` exists before the test runs, it will be deleted
#     and created again.  Therefore, if one wants to preserve the contents of
#     this directory between test runs then one needs to copy the files it
#     contains somewhere else.  This is a good option to use if the commands
#     create intermediate files and one wants to make sure they get deleted
#     before the test cases are run again.  This is also a very useful option
#     to use if multiple tests are defined in the same ``CMakeLists.txt`` file
#     that read/write files with the same name.
#
#   ``FAIL_FAST``
#
#     If specified, then the remaining test commands will be aborted when any
#     test command fails.  Otherwise, all of the test cases will be run.
#
#   ``RUN_SERIAL``
#
#     If specified then no other tests will be allowed to run while this test
#     is running.  This is useful for devices (like CUDA cards) that require
#     exclusive access for processes/threads.  This just sets the CTest test
#     property ``RUN_SERIAL`` using the built-in CMake function
#     ``SET_TESTS_PROPERTIES()``.
#
#   ``COMM [serial] [mpi]``
#
#     If specified, selects if the test will be added in serial and/or MPI
#     mode.  See the ``COMM`` argument in the script
#     `TRIBITS_ADD_TEST()`_ for more details.
#
#   ``OVERALL_NUM_MPI_PROCS <overallNumProcs>``
#
#     If specified, gives the default number of MPI processes that each
#     executable command runs on.  If ``<numProcs>`` is greater than
#     ``${MPI_EXEC_MAX_NUMPROCS}`` then the test will be excluded.  If not
#     specified, then the default number of processes for an MPI build will be
#     ``${MPI_EXEC_DEFAULT_NUMPROCS}``.  For serial builds, this argument is
#     ignored.  This also results in the test property ``PROCESSORS`` being
#     set to ``<overallNumProcs>`` (see `Running multiple tests at the same
#     time (TRIBITS_ADD_ADVANCED_TEST())`_).  **WARNING!** If just running a
#     serial script or other command, then the property ``PROCESSORS`` will
#     still get set to ``${MPI_EXEC_DEFAULT_NUMPROCS}`` so in order to avoid
#     CTest unnecessarily resolving ``${MPI_EXEC_DEFAULT_NUMPROCS}`` processes
#     for a serial non-MPI test, then one must explicitly pass in
#     ``MPI_EXEC_DEFAULT_NUMPROCS 1``!
#
#   ``CATEGORIES <category0> <category1> ...``
#
#     Gives the `Test Test Categories`_ for which this test will be added.
#     See `TRIBITS_ADD_TEST()`_ for more details.
#
#   ``HOST <host0> <host1> ...``
#
#     The list of hosts for which to enable the test (see
#     `TRIBITS_ADD_TEST()`_).
#
#   ``XHOST <host0> <host1> ...``
#
#     The list of hosts for which **not** to enable the test (see
#     `TRIBITS_ADD_TEST()`_).
#
#   ``HOSTTYPE <hosttype0> <hosttype1> ...``
#
#     The list of host types for which to enable the test (see
#     `TRIBITS_ADD_TEST()`_).
#
#   ``XHOSTTYPE <hosttype0> <hosttype1> ...``
#
#     The list of host types for which **not** to enable the test (see
#     `TRIBITS_ADD_TEST()`_).
#
#   ``ENVIRONMENT <var1>=<value1> <var2>=<value2> ..``.
#
#     If passed in, the listed environment variables will be set before
#     calling the test.  This is set using the built-in CTest test property
#     ``ENVIRONMENT``.
#
#   ``TIMEOUT <maxSeconds>``
#
#     If passed in, gives maximum number of seconds the test will be allowed
#     to run before being timed-out (see `TRIBITS_ADD_TEST()`_).  This is for
#     the full CTest test, not individual ``TEST_<idx>`` commands!
#
# .. _TEST_<idx> Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST()):
# 
# **TEST_<idx> Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST())**
#
# Each test command block ``TEST_<idx>`` runs either a package-built test
# executable or some general command executable and is defined as either
# ``EXEC <exeRootName>`` or ``CMND <cmndExec>`` with the arguments:
#
#   ``EXEC <exeRootName> [NOEXEPREFIX] [NOEXESUFFIX] [ADD_DIR_TO_NAME]
#   [DIRECTORY <dir>]``
#
#     If ``EXEC`` is specified, then ``<exeRootName>`` gives the root name of
#     an executable target that will be run as the command.  The full
#     executable name and path is determined in exactly the same way it is in
#     the `TRIBITS_ADD_TEST()`_ function (see `Determining the Executable or
#     Command to Run (TRIBITS_ADD_TEST())`_).  If this is an MPI build, then
#     the executable will be run with MPI using ``NUM_MPI_PROCS <numProcs>``
#     or ``OVERALL_NUM_MPI_PROCS <overallNumProcs>`` (if ``NUM_MPI_PROCS`` is
#     not set for this test case).  If the maximum number of MPI processes
#     allowed is less than this number of MPI processes, then the test will
#     *not* be run.  Note that ``EXEC <exeRootName>`` when ``NOEXEPREFIX`` and
#     ``NOEXESUFFIX`` are specified is basically equivalent to ``CMND
#     <cmndExec>`` except that in an MPI build, ``<exeRootName>`` is always
#     run using MPI.  In this case, one can pass in ``<exeRootName>`` to any
#     command one would like and it will get run with MPI in MPI mode just
#     link any other MPI-enabled built executable.
#
#   ``CMND <cmndExec>``
#
#     If ``CMND`` is specified, then ``<cmndExec>`` gives the executable for a
#     command to be run.  In this case, MPI will never be used to run the
#     executable even when configured in MPI mode
#     (i.e. ``TPL_ENABLE_MPI=ON``).  If one wants to run an arbitrary command
#     using MPI, use ``EXEC <fullPathToCmndExec> NOEXEPREFIX NOEXESUFFIX``
#     instead.
#
# By default, the output (stdout/stderr) for each test command is captured and
# is then echoed to stdout for the overall test.  This is done in order to be
# able to grep the result to determine pass/fail.
#
# Other miscellaneous arguments for each ``TEST_<idx>`` block include:
#
#   ``DIRECTORY <dir>``
#
#     If specified, then the executable is assumed to be in the directory
#     given by relative ``<dir>``.  See `TRIBITS_ADD_TEST()`_.
#
#   ``MESSAGE "<message>"``
#
#     If specified, then the string in ``"<message>"`` will be printed before
#     this test command is run.  This allows adding some documentation about
#     each individual test invocation to make the test output more
#     understandable.
#
#   ``WORKING_DIRECTORY <workingDir>``
#
#     If specified, then the working directory ``<workingDir>`` will be
#     created and the test will be run from within this directory.  If the
#     value ``<workingDir> = TEST_NAME`` is given, then the working directory
#     will be given the name ``${PACKAGE_NAME}_<testName>``.  If the directory
#     ``<workingDir>`` exists before the test runs, it will be deleted and
#     created again.  Therefore, if one wants to preserve the contents of this
#     directory between test runs then one needs to copy the given file
#     somewhere else.  Using a different ``WORKING_DIRECTORY`` for individual
#     test commands allows creating independent working directories for each
#     test case.  This would be useful if a single
#     ``OVERALL_WORKING_DIRECTORY`` was not sufficient for some reason.
#
#   ``NUM_MPI_PROCS <numProcs>``
#
#     If specified, then ``<numProcs>`` is the number of processors used for
#     MPI executables.  If not specified, this will default to
#     ``<overallNumProcs>`` from ``OVERALL_NUM_MPI_PROCS <overallNumProcs>``.
#
#   ``OUTPUT_FILE <outputFile>``
#
#     If specified, then stdout and stderr for the test case will be sent to
#     ``<outputFile>``.  By default, the contents of this file will **also**
#     be printed to STDOUT unless ``NO_ECHO_OUT`` is passed as well.
#
#   ``NO_ECHO_OUTPUT``
#
#     If specified, then the output for the test command will not be echoed to
#     the output for the entire test command.
#
# By default, an atomic test line is assumed to pass if the executable or
# commands returns a non-zero value to the shell.  However, a test case can
# also be defined to pass based on:
#
#   ``PASS_ANY``
#
#     If specified, the test command will be assumed to pass regardless of
#     the return value or any other output.  This would be used when a command
#     that is to follow will determine pass or fail based on output from this
#     command in some way.
#
#   ``PASS_REGULAR_EXPRESSION "<regex>"``
#
#     If specified, the test command will be assumed to pass if it matches the
#     given regular expression.  Otherwise, it is assumed to fail.
#
#   ``PASS_REGULAR_EXPRESSION_ALL "<regex1>" "<regex2>" ... "<regexn>"``
#
#     If specified, the test command will be assumed to pass if the output
#     matches all of the provided regular expressions.  Note that this is not
#     a capability of raw ctest and represents an extension provided by
#     TriBITS.
#
#   ``FAIL_REGULAR_EXPRESSION "<regex>"``
#
#     If specified, the test command will be assumed to fail if it matches the
#     given regular expression.  Otherwise, it is assumed to pass.
#
#   ``STANDARD_PASS_OUTPUT``
#
#     If specified, the test command will be assumed to pass if the string
#     expression "Final Result: PASSED" is found in the output for the test.
#
# All of the arguments for a test block ``TEST_<idx>`` must appear directly
# below their ``TEST_<idx>`` argument and before the next test block (see
# `Argument Parsing and Ordering (TRIBITS_ADD_ADVANCED_TEST())`_).
#
# .. _Overall Pass/Fail (TRIBITS_ADD_ADVANCED_TEST()):
# 
# **Overall Pass/Fail (TRIBITS_ADD_ADVANCED_TEST())**
#
# By default, the overall test will be assumed to pass if it prints::
#
#   "OVERALL FINAL RESULT: TEST PASSED"
#
# However, this can be changed by setting one of the following optional arguments:
#
#   ``FINAL_PASS_REGULAR_EXPRESSION <regex>``
#
#     If specified, the test will be assumed to pass if the output matches
#     ``<regex>``.  Otherwise, it will be assumed to fail.
#
#   ``FINAL_FAIL_REGULAR_EXPRESSION <regex>``
#
#     If specified, the test will be assumed to fail if the output matches
#     ``<regex>``.  Otherwise, it will be assumed to fail.
#
# .. _Argument Parsing and Ordering (TRIBITS_ADD_ADVANCED_TEST()):
# 
# **Argument Parsing and Ordering (TRIBITS_ADD_ADVANCED_TEST())**
#
# The basic tool used for parsing the arguments to this function is the macro
# `PARSE_ARGUMENTS()`_ which has a certain set of behaviors.  The parsing
# using `PARSE_ARGUMENTS()`_ is actually done in two phases.  There is a
# top-level parsing of the "overall" arguments listed in `Overall Arguments
# (TRIBITS_ADD_ADVANCED_TEST())`_ that also pulls out the test blocks.  Then
# there is a second level of parsing using ``PARSE_ARGUMENTS()`` for each of
# the ``TEST_<idx>`` blocks.  Because of this usage, there are a few
# restrictions that one needs to be aware of when using
# ``TRIBITS_ADD_ADVANCED_TEST()``.  This short sections tries to explain the
# behaviors and what is allowed and what is not allowed.
#
# For the most part, the "overall" arguments and the arguments inside of any
# individual ``TEST_<idx>`` blocks can be listed can appear in any order but
# there are restrictions related to the grouping of overall arguments and
# ``TEST_<idx>`` blocks which are as follows:
#
# * The ``<testName>`` argument must be the first listed (it is the only
#   positional argument).
#
# * The test cases ``TEST_<idx>`` must be listed in order (i.e. ``TEST_0
#   ... TEST_1 ...``) and the test cases must be consecutive integers
#   (e.g. can't jump from ``TEST_5`` to ``TEST_7``).
#
# * All of the arguments for a test case must appear directly below its
#   ``TEST_<idx>`` keyword and before the next ``TEST_<idx+1>`` keyword or
#   before any trailing overall keyword arguments.
#
# * None of the overall arguments (e.g. ``CATEGORIES``) can be listed inside
#   of a ``TEST_<idx>`` block but otherwise can be listed before or after all
#   of the ``TEST_<idx>`` blocks.  (NOTE: The current implementation will
#   actually allow overall arguments to be listed after all of the local
#   arguments before the next TEST_<idx> block but this is confusing and will
#   not be allowed in a future implementation).
#
# Other than that, the keyword arguments and options can appear in any order.
#
# .. ToDo: Add some examples of bad argument ordering and what will happen.
#
# .. _Implementation Details (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Implementation Details (TRIBITS_ADD_ADVANCED_TEST())**
#
# Since raw CTest does not support the features provided by this function, the
# way an advanced test is implemented is that a ``cmake -P`` script with the
# name ``${PACKAGE_NAME}_<testName>.cmake`` gets created in the current binary
# directory that then gets added to CTest using::
#
#   ADD_TEST(${PACKAGE_NAME}_<testName>
#     cmake [other options] -P ${PACKAGE_NAME}_<testName>.cmake)
#
# This ``cmake -P`` script then runs the various test cases and checks the
# pass/fail for each case to determine overall pass/fail and implement other
# functionality described above.
#
# .. _Setting Additional Test Properties (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Setting Additional Test Properties (TRIBITS_ADD_ADVANCED_TEST())**
#
# After this function returns, if the test gets added using ``ADD_TEST()``,
# then additional properties can be set and changed using
# ``SET_TESTS_PROPERTIES(${PACKAGE_NAME}_<testName> ...)``.  Therefore, any
# tests properties that are not directly supported by this function and passed
# through the argument list to this wrapper function can be set in the outer
# ``CMakeLists.txt`` file after the call to ``TRIBITS_ADD_ADVANCED_TEST()``.
#
# .. _Running multiple tests at the same time (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Runnning multiple tests at the same time (TRIBITS_ADD_ADVANCED_TEST())**
#
# Just as with `TRIBITS_ADD_TEST()`_, setting ``NUM_MPI_PROCS <numProcs>`` or
# ``OVERALL_NUM_MPI_PROCS <numOverallProcs>`` will set the ``PROCESSORS``
# CTest property to allow CTest to schedule and run mutiple tests at the same
# time when ``'ctest -j<N>'`` is used (see `Running multiple tests at the same
# time (TRIBITS_ADD_TEST())`_).
#
# .. _Disabling Tests Externally (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Disabling Tests Externally (TRIBITS_ADD_ADVANCED_TEST())**
#
# The test can be disabled externally by setting the CMake cache variable
# ``${FULL_TEST_NAME}_DISABLE=TRUE``.  This allows tests to be disabled on a
# case-by-case basis.  The name ``${FULL_TEST_NAME}`` must be the *exact* name
# that shows up in ``ctest -N`` when running the test.
#
# .. _Debugging and Examining Test Generation (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Debugging and Examining Test Generation (TRIBITS_ADD_ADVANCED_TEST())**
#
# In order to see if the test gets added and to debug some issues in test
# creation, one can set the cache variable
# ``${PROJECT_NAME}_VERBOSE_CONFIGURE=ON``.  This will result in the printout
# of some information about the test getting added or not.
#
# Likely the best way to debugging test generation using this function is to
# examine the generated file ``${PACKAGE_NAME}_<testName>.cmake`` in the
# current binary directory (see `Implementation Details
# (TRIBITS_ADD_ADVANCED_TEST())`_) and the generated ``CTestTestfile.cmake``
# file that should list this test case.
#
FUNCTION(TRIBITS_ADD_ADVANCED_TEST TEST_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_ADD_ADVANCED_TEST: ${TEST_NAME_IN}\n")
  ENDIF()

  GLOBAL_SET(TRIBITS_SET_TEST_PROPERTIES_INPUT)

  IF (PACKAGE_NAME)
    SET(TEST_NAME ${PACKAGE_NAME}_${TEST_NAME_IN})
  ELSE()
    SET(TEST_NAME ${TEST_NAME_IN})
  ENDIF()

  #
  # A) Parse the overall arguments and figure out how many tests
  # comands we will have
  #

  # Allow for a maximum of 20 (0 through 19) test commands
  SET(MAX_NUM_TEST_CMND_IDX ${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX})

  SET(TEST_IDX_LIST "")
  FOREACH( TEST_CMND_IDX RANGE ${MAX_NUM_TEST_CMND_IDX})
    LIST( APPEND TEST_IDX_LIST TEST_${TEST_CMND_IDX} )
  ENDFOREACH()
  #PRINT_VAR(TEST_IDX_LIST)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "${TEST_IDX_LIST};OVERALL_WORKING_DIRECTORY;KEYWORDS;COMM;OVERALL_NUM_MPI_PROCS;FINAL_PASS_REGULAR_EXPRESSION;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;FINAL_FAIL_REGULAR_EXPRESSION;TIMEOUT;ENVIRONMENT"
     #options
     "FAIL_FAST;RUN_SERIAL"
     ${ARGN}
     )
  # NOTE: The TIMEOUT argument is not documented on purpose.  I don't want to
  # advertise it!
  
  #
  # B) Add or don't add tests based on a number of criteria
  #

  SET(ADD_THE_TEST FALSE)
  TRIBITS_ADD_TEST_PROCESS_CATEGORIES(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  SET(ADD_THE_TEST FALSE)
  TRIBITS_ADD_TEST_PROCESS_HOST_HOSTTYPE(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  TRIBITS_ADD_TEST_QUERY_DISABLE(DISABLE_THIS_TEST ${TEST_NAME})
  IF (DISABLE_THIS_TEST)
    RETURN()
  ENDIF()

  #
  # C) Determine if we will add the serial or MPI tests based on input COMM
  # and TPL_ENABLE_MPI
  #

  TRIBITS_PROCESS_COMM_ARGS(ADD_SERIAL_TEST  ADD_MPI_TEST  ${PARSE_COMM})
  IF (NOT ADD_SERIAL_TEST AND NOT ADD_MPI_TEST)
    RETURN()
  ENDIF()

  #
  # D) Build the test script
  #

  SET(ADD_THE_TEST TRUE)

  SET(TEST_SCRIPT_STR "")

  APPEND_STRING_VAR( TEST_SCRIPT_STR
    "\n"
    "#\n"
    "# This is a CMake script and must be run as \"cmake -P <SCRIPT_NAME>\"\n"
    "#\n"
    "# NOTE: To see what commands this script runs, run it as:\n"
    "#\n"
    "#    $ cmake -DSHOW_COMMANDS_ONLY=ON -P <SCRIPT_NAME>\n"
    "#\n"
    "\n"
    "#\n"
    "# Variables\n"
    "#\n"
    "\n"
    "SET( TEST_NAME ${TEST_NAME} )\n"
    )

  # Loop through each test case

  SET(NUM_CMNDS 0)
  SET(TEST_EXE_LIST "")

  FOREACH( TEST_CMND_IDX RANGE ${MAX_NUM_TEST_CMND_IDX} )

    IF (NOT PARSE_TEST_${TEST_CMND_IDX} )
      BREAK()
    ENDIF()

    MATH( EXPR NUM_CMNDS ${NUM_CMNDS}+1 )

    # Parse the test command case

    #PRINT_VAR(PARSE_TEST_${TEST_CMND_IDX})

    PARSE_ARGUMENTS(
       #prefix
       PARSE
       #lists
       "EXEC;CMND;ARGS;DIRECTORY;MESSAGE;WORKING_DIRECTORY;OUTPUT_FILE;NUM_MPI_PROCS;PASS_REGULAR_EXPRESSION_ALL;FAIL_REGULAR_EXPRESSION;PASS_REGULAR_EXPRESSION"
       #options
       "NOEXEPREFIX;NOEXESUFFIX;NO_ECHO_OUTPUT;PASS_ANY;STANDARD_PASS_OUTPUT;ADD_DIR_TO_NAME"
       ${PARSE_TEST_${TEST_CMND_IDX}}
       )

    # Write the command

    SET(ARGS_STR ${PARSE_ARGS})
    #PRINT_VAR(ARGS_STR)
    #IF (PARSE_ARGS)
    #  TRIBITS_JOIN_EXEC_PROCESS_SET_ARGS( ARGS_STR ${PARSE_ARGS} )
    #ENDIF()

    IF (PARSE_EXEC)

      LIST( LENGTH PARSE_EXEC PARSE_EXEC_LEN )
      IF (NOT PARSE_EXEC_LEN EQUAL 1)
        MESSAGE(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} EXEC = '${PARSE_EXEC}'"
          " must be a single name.  To add arguments use ARGS <arg1> <arg2> ...." )
      ENDIF()

      TRIBITS_ADD_TEST_GET_EXE_BINARY_NAME( "${PARSE_EXEC}"
        ${PARSE_NOEXEPREFIX} ${PARSE_NOEXESUFFIX} ${PARSE_ADD_DIR_TO_NAME} EXE_BINARY_NAME )

      TRIBITS_ADD_TEST_ADJUST_DIRECTORY( ${EXE_BINARY_NAME} "${PARSE_DIRECTORY}"
        EXECUTABLE_PATH)

      IF (NOT PARSE_NUM_MPI_PROCS)
        SET(PARSE_NUM_MPI_PROCS ${PARSE_OVERALL_NUM_MPI_PROCS})
      ENDIF()
      TRIBITS_ADD_TEST_GET_NUM_PROCS_USED("${PARSE_NUM_MPI_PROCS}" NUM_PROCS_USED)
      IF (NUM_PROCS_USED LESS 0)
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE(STATUS "Skipping adding test because NUM_MPI_PROCS ${PARSE_NUM_MPI_PROCS}"
            " is out of range.")
        ENDIF()
        SET(ADD_THE_TEST FALSE)
      ENDIF()

      IF(ADD_THE_TEST)
        LIST(APPEND TEST_EXE_LIST ${EXECUTABLE_PATH})
      ENDIF()

      TRIBITS_ADD_TEST_GET_TEST_CMND_ARRAY( TEST_CMND_ARRAY
        "${EXECUTABLE_PATH}" "${NUM_PROCS_USED}" ${ARGS_STR} )
      #PRINT_VAR(TEST_CMND_ARRAY)

    ELSEIF (PARSE_CMND)

      LIST( LENGTH PARSE_CMND PARSE_CMND_LEN )
      IF (NOT PARSE_CMND_LEN EQUAL 1)
        MESSAGE(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} CMND = '${PARSE_CMND}'"
          " must be a single command.  To add arguments use ARGS <arg1> <arg2> ...." )
      ENDIF()

      #This allows us to check if a test is requesting more than MPI_EXEC_MAX_NUMPROCS
      TRIBITS_ADD_TEST_GET_NUM_PROCS_USED("${PARSE_OVERALL_NUM_MPI_PROCS}" NUM_PROCS_USED)

      IF (NUM_PROCS_USED LESS 0)
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE(STATUS "Skipping adding test because OVERALL_NUM_MPI_PROCS ${PARSE_OVERALL_NUM_MPI_PROCS}"
            " is out of range.")
        ENDIF()
        SET(ADD_THE_TEST FALSE)
      ENDIF()

      IF(ADD_THE_TEST)
        IF (NOT TRIBITS_ADD_TEST_ADD_TEST_UNITTEST)
          FIND_PROGRAM(CMND_PATH ${PARSE_CMND})
        ELSE()
          SET(CMND_PATH ${PARSE_CMND})
        ENDIF()
        LIST(APPEND TEST_EXE_LIST ${CMND_PATH})
      ENDIF()

      SET( TEST_CMND_ARRAY ${PARSE_CMND} ${ARGS_STR} )
    
    ELSE()

      MESSAGE( FATAL_ERROR
        "Must have EXEC or CMND for TEST_${TEST_CMND_IDX}" )

    ENDIF()

    TRIBITS_JOIN_EXEC_PROCESS_SET_ARGS( TEST_CMND_STR "${TEST_CMND_ARRAY}" )
    #PRINT_VAR(TEST_CMND_STR)

    APPEND_STRING_VAR( TEST_SCRIPT_STR
      "\n"
      "SET( TEST_${TEST_CMND_IDX}_CMND ${TEST_CMND_STR} )\n"
      )
    IF (TRIBITS_ADD_ADVANCED_TEST_UNITTEST)
      GLOBAL_SET(TRIBITS_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX}
        "${TEST_CMND_STR}" )
    ENDIF()

    IF (PARSE_MESSAGE)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_MESSAGE \"${PARSE_MESSAGE}\" )\n"
        )
    ENDIF()

    IF (PARSE_WORKING_DIRECTORY)
      IF ("${PARSE_WORKING_DIRECTORY}" STREQUAL "TEST_NAME")
        SET(PARSE_WORKING_DIRECTORY ${TEST_NAME})
      ENDIF()
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_WORKING_DIRECTORY \"${PARSE_WORKING_DIRECTORY}\" )\n"
        )
    ENDIF()

    IF (PARSE_OUTPUT_FILE)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_OUTPUT_FILE \"${PARSE_OUTPUT_FILE}\" )\n"
        )
    ENDIF()

    IF (PARSE_NO_ECHO_OUTPUT)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_NO_ECHO_OUTPUT \"${PARSE_NO_ECHO_OUTPUT}\" )\n"
        )
    ENDIF()

    # Set up pass/fail

    IF (PARSE_PASS_ANY)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_PASS_ANY TRUE )\n"
        )
    ELSEIF (PARSE_STANDARD_PASS_OUTPUT)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_PASS_REGULAR_EXPRESSION \"End Result: TEST PASSED\" )\n"
        )
    ELSEIF (PARSE_PASS_REGULAR_EXPRESSION)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_PASS_REGULAR_EXPRESSION \"${PARSE_PASS_REGULAR_EXPRESSION}\" )\n"
        )
    ELSEIF (PARSE_PASS_REGULAR_EXPRESSION_ALL)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_PASS_REGULAR_EXPRESSION_ALL "
        )
      FOREACH(REGEX_STR ${PARSE_PASS_REGULAR_EXPRESSION_ALL})
        APPEND_STRING_VAR( TEST_SCRIPT_STR
          "\"${REGEX_STR}\" "
          )
      ENDFOREACH()
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        ")\n"
        )
    ENDIF()

  ENDFOREACH()

  # ToDo: Verify that TEST_${MAX_NUM_TEST_CMND_IDX}+1 does *not* exist!

  IF (PARSE_OVERALL_WORKING_DIRECTORY)
    IF ("${PARSE_OVERALL_WORKING_DIRECTORY}" STREQUAL "TEST_NAME")
      SET(PARSE_OVERALL_WORKING_DIRECTORY ${TEST_NAME})
    ENDIF()
  ENDIF()

  APPEND_STRING_VAR( TEST_SCRIPT_STR
    "\n"
    "SET(NUM_CMNDS ${NUM_CMNDS})\n"
    "\n"
    "SET(OVERALL_WORKING_DIRECTORY \"${PARSE_OVERALL_WORKING_DIRECTORY}\")\n"
    "\n"
    "SET(FAIL_FAST ${PARSE_FAIL_FAST})\n"
    "\n"
    "#\n"
    "# Test invocation\n"
    "#\n"
    "\n"
    "SET(CMAKE_MODULE_PATH ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_UTILS_DIR})\n"
    "\n"
    "INCLUDE(DriveAdvancedTest)\n"
    "\n"
    "DRIVE_ADVANCED_TEST()\n"
    )

  IF (TRIBITS_ADD_ADVANCED_TEST_UNITTEST)
    GLOBAL_SET(TRIBITS_ADD_ADVANCED_TEST_NUM_CMNDS ${NUM_CMNDS})
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(TEST_SCRIPT_STR)
  ENDIF()

  # Write the script file

  SET(TEST_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.cmake")

  IF (ADD_THE_TEST AND NOT TRIBITS_ADD_ADVANCED_TEST_SKIP_SCRIPT)

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("\nWriting file \"${TEST_SCRIPT_FILE}\" ...")
    ENDIF()
  
    FILE( WRITE "${TEST_SCRIPT_FILE}"
      "${TEST_SCRIPT_STR}" )

  ENDIF()

  #
  # F) Set the CTest test to run the new script
  #

  IF (ADD_THE_TEST)

    IF(NOT TRIBITS_ADD_TEST_ADD_TEST_UNITTEST)
      # Tell CTest to run our script for this test.  Pass the test-type
      # configuration name to the script in the TEST_CONFIG variable.
      ADD_TEST( ${TEST_NAME}
        ${CMAKE_COMMAND} "-DTEST_CONFIG=\${CTEST_CONFIGURATION_TYPE}"
        -P "${TEST_SCRIPT_FILE}")
    ENDIF()

    LIST(REMOVE_DUPLICATES TEST_EXE_LIST)
    TRIBITS_SET_TEST_PROPERTY(${TEST_NAME} PROPERTY REQUIRED_FILES ${TEST_EXE_LIST})
  
    IF(PARSE_RUN_SERIAL)
      TRIBITS_SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES RUN_SERIAL ON)
    ENDIF()
    
    TRIBITS_PRIVATE_ADD_TEST_ADD_LABEL_AND_KEYWORDS(${TEST_NAME})
    
    #This if clause will set the number of PROCESSORS to reserve during testing
    #to the number requested for the test.
    IF(NUM_PROCS_USED)
      TRIBITS_SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES
        PROCESSORS "${NUM_PROCS_USED}")
    ENDIF()
    
    IF (PARSE_FINAL_PASS_REGULAR_EXPRESSION)
      TRIBITS_SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        PASS_REGULAR_EXPRESSION "${PARSE_FINAL_PASS_REGULAR_EXPRESSION}" )
    ELSEIF (PARSE_FINAL_FAIL_REGULAR_EXPRESSION)
      TRIBITS_SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        FAIL_REGULAR_EXPRESSION "${PARSE_FINAL_FAIL_REGULAR_EXPRESSION}" )
    ELSE()
      TRIBITS_SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        PASS_REGULAR_EXPRESSION "OVERALL FINAL RESULT: TEST PASSED" )
    ENDIF()
    
    IF (PARSE_TIMEOUT)
      TRIBITS_SCALE_TIMEOUT(${PARSE_TIMEOUT} TIMEOUT_USED)
      TRIBITS_SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES TIMEOUT ${TIMEOUT_USED})
    ENDIF()
  
    IF (PARSE_ENVIRONMENT)
      TRIBITS_SET_TEST_PROPERTY(${TEST_NAME} PROPERTY ENVIRONMENT ${PARSE_ENVIRONMENT})
    ENDIF()

  ENDIF()

ENDFUNCTION()

# PERFORMANCE NOTES:
#
# We might be able to improve the performance of the parsing by limiting the
# number of TEST_<I> blocks up front by setting a varible that will fix it.
# This might just be set as a local variable in the CMakeLists.txt file where
# this function is called from.  The other option is to just read through the
# input arguments to parse first and look for the highest TEST_<idx> and use
# that instead to build the list of tests.  This would just be a linear search
# it could be a big pay off for very long lists.
