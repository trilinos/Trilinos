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

include(TribitsCMakePolicies  NO_POLICY_SCOPE)
include(TribitsAddTestHelpers)


# @FUNCTION: tribits_add_test()
#
# Add a test or a set of tests for a single executable or command using CTest
# ``add_test()``.
#
# Usage::
#
#   tribits_add_test(
#     <exeRootName>  [NOEXEPREFIX]  [NOEXESUFFIX]
#     [NAME <testName> | NAME_POSTFIX <testNamePostfix>]
#     [DIRECTORY <directory>]
#     [ADD_DIR_TO_NAME]
#     [RUN_SERIAL]
#     [ARGS "<arg0> <arg1> ..." "<arg2> <arg3> ..." ...
#       | POSTFIX_AND_ARGS_0 <postfix0> <arg0> <arg1> ...
#         POSTFIX_AND_ARGS_1 ... ]
#     [COMM [serial] [mpi]]
#     [NUM_MPI_PROCS <numMpiProcs>]
#     [NUM_TOTAL_CORES_USED <numTotalCoresUsed>]
#     [CATEGORIES <category0>  <category1> ...]
#     [HOST <host0> <host1> ...]
#     [XHOST <host0> <host1> ...]
#     [HOSTTYPE <hosttype0> <hosttype1> ...]
#     [XHOSTTYPE <hosttype0> <hosttype1> ...]
#     [EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...]
#     [DISABLED <messageWhyDisabled>]
#     [STANDARD_PASS_OUTPUT
#       | PASS_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...]
#     [FAIL_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...]
#     [WILL_FAIL]
#     [ENVIRONMENT <var0>=<value0> <var1>=<value1> ...]
#     [TIMEOUT <maxSeconds>]
#     [LIST_SEPARATOR <sep>]
#     [ADDED_TESTS_NAMES_OUT <testsNames>]
#     )
#
# The tests are only added if tests are enabled for the package
# (i.e. `${PACKAGE_NAME}_ENABLE_TESTS`_ ``= ON``).  (NOTE: A more efficient
# way to optionally enable tests or examples is to put them in a ``test/`` or
# ``example/`` subdir and then include that subdir with
# `tribits_add_test_directories()`_ or `tribits_add_example_directories()`_,
# respectively.)
#
# *Sections:*
#
# * `Formal Arguments (tribits_add_test())`_
# * `Determining the Executable or Command to Run (tribits_add_test())`_
# * `Determining the Full Test Name (tribits_add_test())`_
# * `Adding Multiple Tests  (tribits_add_test())`_
# * `Determining Pass/Fail (tribits_add_test())`_
# * `Setting additional test properties (tribits_add_test())`_
# * `Running multiple tests at the same time (tribits_add_test())`_
# * `Setting timeouts for tests (tribits_add_test())`_
# * `Debugging and Examining Test Generation (tribits_add_test())`_
# * `Disabling Tests Externally (tribits_add_test())`_
# * `Adding extra commandline arguments externally (tribits_add_test())`_
#
# .. _Formal Arguments (tribits_add_test()):
#
# **Formal Arguments (tribits_add_test())**
#
#   ``<exeRootName>``
#
#     The name of the executable or path to the executable to run for the test
#     (see `Determining the Executable or Command to Run
#     (tribits_add_test())`_).  This name is also the default root name for
#     the test (see `Determining the Full Test Name (tribits_add_test())`_).
#
#   ``NOEXEPREFIX``
#
#    If specified, then the prefix ``${PACKAGE_NAME}_`` is assumed **not** to
#    be prepended to ``<exeRootName>`` (see `Determining the Executable or
#    Command to Run (tribits_add_test())`_).
#
#   ``NOEXESUFFIX``
#
#      If specified, then the postfix
#      ``${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}`` is assumed **not** to be
#      post-pended to ``<exeRootName>`` (except on Windows platforms, see
#      `Determining the Executable or Command to Run (tribits_add_test())`_).
#
#   ``NAME <testRootName>``
#
#     If specified, gives the root name of the test.  If not specified, then
#     ``<testRootName>`` is taken to be ``<exeRootName>``.  The actual test
#     name passed to ``add_test()`` will always be prefixed as
#     ``${PACKAGE_NAME}_<testRootName>``.  The main purpose of this argument
#     is to allow multiple tests to be defined for the same executable.  CTest
#     requires all test names to be globally unique in a single project.  See
#     `Determining the Full Test Name (tribits_add_test())`_.
#
#   ``NAME_POSTFIX <testNamePostfix>``
#
#     If specified, gives a postfix that will be added to the standard test
#     name based on ``<exeRootName>`` (appended as ``_<NAME_POSTFIX>``).  If
#     the ``NAME <testRootName>`` argument is given, this argument is ignored.
#     See `Determining the Full Test Name (tribits_add_test())`_.
#
#   ``DIRECTORY <dir>``
#
#     If specified, then the executable is assumed to be in the directory
#     given by ``<dir>``.  The directory ``<dir>`` can either be a relative or
#     absolute path.  If not specified, the executable is assumed to be in the
#     current binary directory ``${CMAKE_CURRENT_BINARY_DIR}``.  See
#     `Determining the Executable or Command to Run (tribits_add_test())`_.
#
#   ``ADD_DIR_TO_NAME``
#
#     If specified, then the directory name that this test resides in will be
#     added into the name of the test after the package name is added and
#     before the root test name (see `Determining the Full Test Name
#     (tribits_add_test())`_).  The directory name will have the package's
#     base directory stripped off so only the unique part of the test
#     directory will be used.  All directory separators ``"/"`` will be
#     changed into underscores ``"_"``.
#
#   ``RUN_SERIAL``
#
#     If specified, then no other tests will be allowed to run while this test
#     is running. This is useful for devices (like CUDA GPUs) that require
#     exclusive access for processes/threads.  This just sets the CTest test
#     property ``RUN_SERIAL`` using the built-in CMake function
#     ``set_tests_properties()``.  Also, the addition of the ``RUN_SERIAL``
#     test property can be triggered by (the user) setting the global cache
#     variable ``<fullTestName>_SET_RUN_SERIAL=ON``.  NOTE: If ``RUN_SERIAL``
#     is passed in but ``<fullTestName>_SET_RUN_SERIAL=OFF`` (or any value
#     evaluating to ``FALSE``), then the ``RUN_SERIAL`` test property will
#     **NOT** be set on the added test(s).
#
#   ``ARGS "<arg0> <arg1> ..." "<arg2> <arg3> ..." ...``
#
#     If specified, then a set of arguments can be passed in quotes.  If
#     multiple groups of arguments are passed in different quoted clusters of
#     arguments then a different test will be added for each set of arguments.
#     In this way, many different tests can be added for a single executable
#     in a single call to this function.  Each of these separate tests will be
#     named ``<fullTestName>_xy`` where ``xy`` = ``00``, ``01``, ``02``, and so
#     on.  **WARNING:** When defining multiple tests it is preferred to use the
#     ``POSTFIX_AND_ARGS_<IDX>`` form instead.  **WARNING:** Multiple
#     arguments passed to a single test invocation must be quoted or multiple
#     tests taking single arguments will be created instead!  See `Adding
#     Multiple Tests (tribits_add_test())`_ for more details and examples.
#
#   ``POSTFIX_AND_ARGS_<IDX> <postfix> <arg0> <arg1> ...``
#
#     If specified, gives a sequence of sets of test postfix names and
#     arguments lists for different tests (up to ``POSTFIX_AND_ARGS_19``).
#     For example, a set of three different tests with argument lists can be
#     specified as::
#
#       POSTIFX_AND_ARGS_0 postfix0 --arg1 --arg2=dummy
#       POSTIFX_AND_ARGS_1 postfix1  --arg2=fly
#       POSTIFX_AND_ARGS_2 postfix2  --arg2=bags
#
#     This will create three different test cases with the postfix names
#     ``postfix0``, ``postfix1``, and ``postfix2``.  The indexes must be
#     consecutive starting a ``0`` and going up to (currently) ``19``.  The
#     main advantages of using these arguments instead of just ``ARGS`` are
#     that one can give a meaningful name to each test case and one can
#     specify multiple arguments without having to quote them and one can
#     allow long argument lists to span multiple lines.  See `Adding Multiple
#     Tests (tribits_add_test())`_ for more details and examples.
#
#     Note that one of the `<postfix>` arguments can be empty, in which case
#     the base test name is not appended so::
#
#       POSTIFX_AND_ARGS_0 "" --arg1 --arg2=dummy
#
#     would create one test without appending the test name.
#
#   ``COMM [serial] [mpi]``
#
#     If specified, determines if the test will be added in serial and/or MPI
#     mode.  If the ``COMM`` argument is missing, the test will be added in
#     both serial and MPI builds of the code.  That is if ``COMM mpi`` is
#     passed in, then the test will **not** be added if
#     ``TPL_ENABLE_MPI=OFF``.  Likewise, if ``COMM serial`` is passed in, then
#     the test will **not** be added if ``TPL_ENABLE_MPI=ON``.  If ``COMM
#     serial mpi`` or ``COMM mpi serial`` is passed in, then the value of
#     ``TPL_ENABLE_MPI`` does not determine if the test is added or not.
#
#   ``NUM_MPI_PROCS <numMpiProcs>``
#
#     If specified, gives the number of MPI processes used to run the test
#     with the MPI exec program ``${MPI_EXEC}``.  If ``<numMpiProcs>`` is
#     greater than ``${MPI_EXEC_MAX_NUMPROCS}`` then the test will be
#     excluded.  If not specified, then the default number of processes for an
#     MPI build (i.e. ``TPL_ENABLE_MPI=ON``) will be
#     ``${MPI_EXEC_DEFAULT_NUMPROCS}``.  For serial builds
#     (i.e. ``TPL_ENABLE_MPI=OFF``), passing in a value ``<numMpiProcs>`` >
#     ``1`` will cause the test to not be added.  The value ``<numMpiProcs>``
#     will also be set as the built-in test property ``PROCESSORS`` if
#     ``NUM_TOTAL_CORES_USED`` is not specified.
#
#   ``NUM_TOTAL_CORES_USED <numTotalCoresUsed>``
#
#     If specified, gives the total number of processes or cores that is
#     reported to CTest as the built-in CTest ``PROCESSORS`` property.  If
#     this is not specified, then ``PROCESSORS`` is specified by the argument
#     ``NUM_MPI_PROCS <numMpiProcs>``.  This argument is used for test
#     scripts/executables that use more cores than MPI processes
#     (i.e. ``<numMpiProcs>``) and its only purpose is to inform CTest and
#     TriBITS of the maximum number of processes or cores that are used by the
#     underlying test executable/script.  When specified, if
#     ``<numTotalCoresUsed>`` is greater than ``${MPI_EXEC_MAX_NUMPROCS}``,
#     then the test will not be added.  Otherwise, the CTest property
#     ``PROCESSORS`` is set to ``<numTotalCoresUsed>`` so that CTest knows how
#     to best schedule the test w.r.t. other tests on a given number of
#     available processes.  See `Running multiple tests at the same time
#     (tribits_add_test())`_.
#
#   ``CATEGORIES <category0> <category1> ...``
#
#     If specified, gives the specific categories of the test.  Valid test
#     categories include ``BASIC``, ``CONTINUOUS``, ``NIGHTLY``, ``HEAVY``
#     and ``PERFORMANCE``.  If not specified, the default category is
#     ``BASIC``.  When the test category does not match
#     ``${PROJECT_NAME}_TEST_CATEGORIES``, then the test is **not** added.
#     When ``CATEGORIES`` contains ``BASIC`` it will match
#     ``${PROJECT_NAME}_TEST_CATEGORIES`` equal to ``CONTINUOUS``,
#     ``NIGHTLY``, and ``HEAVY``.  When ``CATEGORIES`` contains
#     ``CONTINUOUS`` it will match ``${PROJECT_NAME}_TEST_CATEGORIES`` equal
#     to ``CONTINUOUS``, ``NIGHTLY``, and ``HEAVY``.  When ``CATEGORIES``
#     contains ``NIGHTLY`` it will match ``${PROJECT_NAME}_TEST_CATEGORIES``
#     equal to ``NIGHTLY`` and ``HEAVY``.  When ``CATEGORIES`` contains
#     ``PERFORMANCE`` it will match
#     ``${PROJECT_NAME}_TEST_CATEGORIES=PERFORMANCE`` only.
#
#   ``HOST <host0> <host1> ...``
#
#     If specified, gives a list of hostnames where the test will be included.
#     The current hostname is determined by the built-in CMake command
#     ``site_name(${PROJECT_NAME}_HOSTNAME)``.  On Linux/Unix systems, this is
#     typically the value returned by ``uname -n``.  If this list is given,
#     the value of ``${${PROJECT_NAME}_HOSTNAME}`` must equal one of the
#     listed host names ``<hosti>`` or test will **not** be added.  The value
#     of ``${PROJECT_NAME}_HOSTNAME`` gets printed out in the TriBITS cmake
#     output under the section ``Probing the environment`` (see `Full
#     Processing of TriBITS Project Files`_).
#
#   ``XHOST <host0> <host1> ...``
#
#     If specified, gives a list of hostnames (see ``HOST`` argument) on which
#     the test will **not** be added.  This check is performed after the check
#     for the hostnames in the ``HOST`` list if it should exist.  Therefore,
#     this exclusion list overrides the ``HOST`` inclusion list.
#
#   ``HOSTTYPE <hosttype0> <hosttype1> ...``
#
#     If specified, gives the names of the host system type (given by the
#     built-in CMake cache variable ``CMAKE_HOST_SYSTEM_NAME`` which is
#     printed in the TriBITS cmake configure output in the section ``Probing
#     the environment``) for which the test is allowed to be added.  If
#     ``HOSTTYPE`` is specified and ``CMAKE_HOST_SYSTEM_NAME`` is not equal to
#     one of the values of ``<hosttypei>``, then the test will **not** be
#     added.  Typical host system type names include ``Linux``, ``Darwin``,
#     ``Windows``, etc.
#
#   ``XHOSTTYPE <hosttype0> <hosttype1> ...``
#
#     If specified, gives the names of the host system type (see the
#     ``HOSTTYPE`` argument above) for which **not** to include the test on.
#     This check is performed after the check for the host system names in the
#     ``HOSTTYPE`` list if it should exist.  Therefore, this exclusion list
#     overrides the ``HOSTTYPE`` inclusion list.
#
#   ``EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...``
#
#     If specified, gives the names of CMake variables that must evaluate to
#     true, or the test will not be added.
#
#   ``DISABLED <messageWhyDisabled>``
#
#     If ``<messageWhyDisabled>`` is non-empty and does not evaluate to FALSE,
#     then the test will be added by ``add_test()`` (so CTest will see it) but
#     the ctest test property ``DISABLED`` will be set.  Therefore, CTest will
#     not run the test and will instead list it as "Not Run" when tests are
#     run locally and when submitting test results to CDash (with test details
#     "Not Run (Disabled)").  Also, the message ``<messageWhyDisabled>`` will
#     be printed to STDOUT by cmake after the line stating the test was added
#     when ``${PROJECT_NAME}_TRACE_ADD_TEST=ON`` is set.  If
#     ``<messageWhyDisabled>`` evaluates to FALSE in CMake (e.g. "FALSE",
#     "false", "NO", "no", "0", "", etc.), then the ``DISABLED`` property will
#     not be set.  This property can also be set with the CMake cache var ``-D
#     <fullTestName>_SET_DISABLED_AND_MSG="<msgSetByVar>"`` and in fact that
#     var will override the value of ``<messageWhyDisabled>`` passed in here
#     (if ``<msgSetByVar>`` is non-empty).  This allows a user to enable tests
#     that are disabled in the CMakeList.txt files using this input.
#
#   ``STANDARD_PASS_OUTPUT``
#
#     If specified, then the standard test output string ``End Result: TEST
#     PASSED`` is grepped in the test stdout for to determine success.  This
#     is needed for MPI tests on some platforms since the return value from
#     MPI executables is unreliable.  This is set using the built-in CTest
#     property ``PASS_REGULAR_EXPRESSION``.
#
#   ``PASS_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...``
#
#     If specified, then the test will be assumed to pass only if one of the
#     regular expressions ``<regex0>``, ``<regex1>`` etc. match the output
#     send to stdout.  Otherwise, the test will fail.  This is set using the
#     built-in CTest property ``PASS_REGULAR_EXPRESSION``.  Consult standard
#     CMake documentation for full behavior.  TIPS: Replace ';' with '[;]' or
#     CMake will interpret this as an array element boundary.  To match '.',
#     use '[.]'.
#
#   ``FAIL_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...``
#
#     If specified, then a test will be assumed to fail if one of the regular
#     expressions ``<regex0>``, ``<regex1>`` etc. match the output send to
#     stdout.  Otherwise, the test will pass.  This is set using the built-in
#     CTest property ``FAIL_REGULAR_EXPRESSION``.  Consult standard CMake
#     documentation for full behavior (and see above tips for
#     ``PASS_REGULAR_EXPRESSION``).
#
#   ``WILL_FAIL``
#
#     If passed in, then the pass/fail criteria will be inverted.  This is set
#     using the built-in CTest property ``WILL_FAIL``.  Consult standard CMake
#     documentation for full behavior.
#
#   ``ENVIRONMENT "<var1>=<value1>" "<var2>=<value2>" ...``.
#
#     If passed in, the listed environment variables will be set by CTest
#     before calling the test.  This is set using the built-in CTest test
#     property ``ENVIRONMENT``.  Note, if the env var values contain
#     semi-colons ``';'``, then replace the semi-colons ``';'`` with another
#     separator ``'<sep>'`` and pass in ``LIST_SEPARATOR <sep>`` so ``<sep>``
#     will be replaced with ``';'`` at point of usage.  If the env var values
#     contain any spaces, also quote the entire variable/value pair as
#     ``"<vari>=<valuei>"``.  For example, the env var and value
#     ``my_env_var="arg1 b;arg2;I have spaces"`` would need to be passed as
#     ``"my_env_var=arg1 b<sep>arg2<sep>I have spaces"``.
#
#   ``TIMEOUT <maxSeconds>``
#
#     If passed in, gives maximum number of seconds the test will be allowed
#     to run before being timed-out and killed.  This sets the CTest property
#     ``TIMEOUT``.  The value ``<maxSeconds>`` will be scaled by the value of
#     `${PROJECT_NAME}_SCALE_TEST_TIMEOUT`_.  See `Setting timeouts for tests
#     (tribits_add_test())`_ for more details.
#
#     **WARNING:** Rather than just increasing the timeout for an expensive
#     test, please try to either make the test run faster or relegate the test
#     to being run less often (i.e. set ``CATEGORIES NIGHTLY`` or even
#     ``HEAVY`` for extremely expensive tests).  Expensive tests are one of
#     the worse forms of technical debt that a project can have!
#
#   ``LIST_SEPARATOR <sep>``
#
#     String used as placeholder for the semi-colon char ``';'`` in order to
#     allow pass-through.  For example, if arguments to the ``ARGS`` or
#     ``ENVIRONMENT`` need to use semi-colons, then replace ``';'`` with
#     ``'<semicolon>'`` (for example) such as with
#     ``"somearg=arg1<semicolon>arg2"``, then at the point of usage,
#     ``'<semicolon>'`` will be replaced with ``';'`` and it will be passed to
#     the final command as ``"somearg=arg1;arg2"`` (with as many preceding
#     escape backslashes ``'\'`` in front of ``';'`` as is needed for the
#     given usage context).
#
#   ``ADDED_TESTS_NAMES_OUT <testsNames>``
#
#     If specified, then on output the variable ``<testsNames>`` will be set
#     with the name(S) of the tests passed to ``add_test()``.  If more than
#     one test is added, then this will be a list of test names.  Having this
#     name allows the calling ``CMakeLists.txt`` file access and set
#     additional test properties (see `Setting additional test properties
#     (tribits_add_test())`_).
#
# In the end, this function just calls the built-in CMake commands
# ``add_test(${TEST_NAME} ...)`` and ``set_tests_properties(${TEST_NAME}
# ...)`` to set up a executable process for ``ctest`` to run, determine
# pass/fail criteria, and set some other test properties.  Therefore, this
# wrapper function does not provide any fundamentally new features that are
# not already available in the basic usage if CMake/CTest.  However, this
# wrapper function takes care of many of the details and boiler-plate CMake
# code that it takes to add such a test (or tests) and enforces consistency
# across a large project for how tests are defined, run, and named (to avoid
# test name clashes).
#
# If more flexibility or control is needed when defining tests, then the
# function `tribits_add_advanced_test()`_ should be used instead.
#
# In the following subsections, more details on how tests are defined and run
# is given.
#
# .. _Determining the Executable or Command to Run (tribits_add_test()):
#
# **Determining the Executable or Command to Run (tribits_add_test())**
#
# This function is primarily designed to make it easy to run tests for
# executables built using the function `tribits_add_executable()`_.  To set up
# tests to run arbitrary executables, see below.
#
# By default, the executable to run is determined by first getting the
# executable name which by default is assumed to be::
#
#  <fullExeName> =
#    ${PACKAGE_NAME}_<exeRootName>${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}
#
# which is (by no coincidence) identical to how it is selected in
# `tribits_add_executable()`_ (see `Executable and Target Name
# (tribits_add_executable())`_).  This name can be altered by passing in
# ``NOEXEPREFIX``, ``NOEXESUFFIX``, and ``ADD_DIR_TO_NAME`` as described in
# `Executable and Target Name (tribits_add_executable())`_.
#
# By default, this executable is assumed to be in the current CMake binary
# directory ``${CMAKE_CURRENT_BINARY_DIR}`` but the directory location can be
# changed using the ``DIRECTORY <dir>`` argument.
#
# If an arbitrary executable is to be run (i.e. not build inside of the
# project), then pass in ``NOEXEPREFIX`` and ``NOEXESUFFIX`` and set
# ``<exeRootName>`` to the relative or absolute path of the executable to be
# run.  If ``<exeRootName>`` is not an absolute path, then
# ``${CMAKE_CURRENT_BINARY_DIR}/<exeRootName>`` is set as the executable to
# run in this case.
#
# NOTE: On native Windows platforms, the ``NOEXESUFFIX`` will still allow
# CTest to run executables that have the ``*.exe`` suffix.
#
# Whatever executable path is specified using this logic, if the executable is
# not found, then when ``ctest`` goes to run the test, it will mark it as
# ``NOT RUN``.
#
# .. _Determining the Full Test Name (tribits_add_test()):
#
# **Determining the Full Test Name (tribits_add_test())**
#
# By default, the base test name is selected to be::
#
#   <fullTestName> = ${PACKAGE_NAME}_<exeRootName>
#
# If ``NAME <testRootName>`` is passed in, then ``<testRootName>`` is used
# instead of ``<exeRootName>`` above.
#
# If ``NAME_POSTFIX <testNamePostfix>`` is passed in, then the base test name
# is selected to be::
#
#   <fullTestName> = ${PACKAGE_NAME}_<exeRootName>_<testNamePostfix>
#
# If ``ADD_DIR_TO_NAME`` is passed in, then the directory name relative to the
# package base directory is added to the name as well to help disambiguate the
# test name (see the above).
#
# Let the test name determined as described above be ``<fullTestName>``.  If
# no arguments or only a single set of arguments are passed in through
# ``ARGS``, then this is the test name actually passed in to ``add_test()``.
# If multiple tests are defined, then this name becomes the base test name for
# each of the tests (see `Adding Multiple Tests (tribits_add_test())`_).
#
# Finally, for any test that gets defined, if MPI is enabled
# (i.e. ``TPL_ENABLE_MPI=ON``), then the terminal suffix
# ``_MPI_${NUM_MPI_PROCS}`` will be added to the end of the test name (even
# for multiple tests).  No such prefix is added for the serial case
# (i.e. ``TPL_ENABLE_MPI=OFF``).
#
# .. _Adding Multiple Tests  (tribits_add_test()):
#
# **Adding Multiple Tests  (tribits_add_test())**
#
# Using this function, one can add executable arguments and can even add
# multiple tests in one of two ways.  One can either pass in one or more
# **quoted** clusters of arguments using::
#
#   ARGS "<arg0> <arg1> ..." "<arg2> <arg3> ..." ...
#
# or can pass in an explicit test name postfix and arguments with::
#
#   POSTFIX_AND_ARGS_0 <postfix0> <arg0> <arg1> ...
#   POSTFIX_AND_ARGS_1 <postfix1> <arg2> ...
#   ...
#
# If only one short set of arguments needs to be passed in, then passing::
#
#   ARGS "<arg0> <arg1>"
#
# may be preferable since it will not add any postfix name to the test.  To
# add more than one test case using ``ARGS``, one will use more than one
# quoted set of arguments such as with::
#
#   ARGS "<arg0> <arg1>" "<arg2> <arg2>"
#
# which creates 2 tests with the names ``<fullTestName>_00`` passing
# arguments ``"<arg0> <arg1>"`` and ``<fullTestName>_01`` passing arguments
# ``"<arg2> <arg3>"``.  However, when passing multiple sets of arguments it is
# preferable to **not** use ``ARGS`` but instead use::
#
#   POSTFIX_AND_ARGS_0 test_a <arg0> <arg1>
#   POSTFIX_AND_ARGS_1 test_b <arg2> <arg2>
#
# which also creates the same 2 tests but now with the improved names
# ``<fullTestName>_test_a`` passing arguments ``"<arg0> <arg1>"`` and
# ``<fullTestName>_test_b`` passing arguments ``"<arg2> <arg3>"``.  In this way,
# the individual tests can be given more understandable names.
#
# The other advantage of the ``POSTFIX_AND_ARGS_<IDX>`` form is that the
# arguments ``<arg0>``, ``<arg1>``, ... do not need to be quoted and can
# therefore be extended over multiple lines like::
#
#   POSTFIX_AND_ARGS_0 long_args --this-is-the-first-long-arg=very
#     --this-is-the-second-long-arg=verylong
#
# If one does not use quotes when using ``ARGS`` one will actually get more
# than one test.  For example, if one passes in::
#
#   ARGS --this-is-the-first-long-arg=very
#     --this-is-the-second-long-arg=verylong
#
# one actually gets two tests, not one test.  This is a common mistake that
# people make when using the ``ARGS`` form of passing arguments.  This can't
# be fixed or it will break backward compatibility.  If this could be designed
# fresh, the ``ARGS`` argument would only create a single test and the
# arguments would not be quoted.
#
# .. _Determining Pass/Fail (tribits_add_test()):
#
# **Determining Pass/Fail (tribits_add_test())**
#
# The only means to determine pass/fail is to use the built-in CTest
# properties ``PASS_REGULAR_EXPRESSION`` and ``FAIL_REGULAR_EXPRESSION`` which
# can only grep the test's STDOUT/STDERR or to check for a 0 return value (or
# invert these using ``WILL_FAIL``).  For simple tests, that is enough.
# However, for more complex executables, one may need to examine one or more
# output files to determine pass/fail.  Raw CMake/CTest cannot do this.  In
# this case, one should use `tribits_add_advanced_test()`_ instead to add the
# test.
#
# .. _Setting additional test properties (tribits_add_test()):
#
# **Setting additional test properties (tribits_add_test())**
#
# After this function returns, any tests that get added using ``add_test()``
# can have additional properties set and changed using
# ``set_tests_properties()``.  Therefore, any tests properties that are not
# directly supported and passed through this wrapper function can be set in
# the outer ``CMakeLists.txt`` file after the call to ``tribits_add_test()``.
#
# If tests are added, then the names of those tests will be returned in the
# variable ``ADDED_TESTS_NAMES_OUT <testsNames>``.  This can be used, for
# example, to override the ``PROCESSORS`` property for the tests with::
#
#   tribits_add_test( someTest ...
#     ADDED_TESTS_NAMES_OUT  someTest_TEST_NAME )
#
#   if (someTest_TEST_NAME)
#     set_tests_properties( ${someTest_TEST_NAME}
#       PROPERTIES ATTACHED_FILES someTest.log )
#   endif()
#
# where the test writes a log file ``someTest.log`` that we want to submit to
# CDash also.
#
# This approach will work no matter what TriBITS names the individual test(s)
# or whether the test(s) are added or not (depending on other arguments like
# ``COMM``, ``XHOST``, etc.).
#
# The following built-in CTest test properties are set through `Formal
# Arguments (tribits_add_test())`_ or are otherwise automatically set by this
# function and should **NOT** be overridden by direct calls to
# ``set_tests_properties()``: ``ENVIRONMENT``, ``FAIL_REGULAR_EXPRESSION``,
# ``LABELS``, ``PASS_REGULAR_EXPRESSION``, ``RUN_SERIAL``, ``TIMEOUT``, and
# ``WILL_FAIL``.
#
# However, generally, other built-in CTest test properties can be set after
# the test is added like show above.  Examples of test properties that can be
# set using direct calls to ``set_tests_properties()`` include
# ``ATTACHED_FILES``, ``ATTACHED_FILES_ON_FAIL``, ``COST``, ``DEPENDS``,
# ``MEASUREMENT``, ``RESOURCE_LOCK`` and ``WORKING_DIRECTORY``.
#
# For example, one can set a dependency between two tests using::
#
#   tribits_add_test( test_a [...]
#      ADDED_TESTS_NAMES_OUT  test_a_TEST_NAME )
#   
#   tribits_add_test( test_b [...]
#      ADDED_TESTS_NAMES_OUT  test_z_TEST_NAME )
#   
#   if (test_a_TEST_NAME AND test_b_TEST_NAME)
#     set_tests_properties(${test_b_TEST_NAME}
#       PROPERTIES DEPENDS ${test_a_TEST_NAME})
#   endif()
#
# This ensures that test ``test_b`` will always be run after ``test_a`` if
# both tests are run by CTest.
#
# .. _Running multiple tests at the same time (tribits_add_test()):
#
# **Running multiple tests at the same time (tribits_add_test())**
#
# By default, CTest will run as many tests defined with ``add_test()`` at same
# time as it can according to its parallel level (e.g. ``'ctest -j<N>'`` or
# the CTest property ``CTEST_PARALLEL_LEVEL``).  For example, when raw
# ``'ctest -j10'`` is run, CTest will run multiple tests at the same time to
# try to make usage of 10 processors/cores.  If all of the defined tests only
# used one process (which is assumed by default except for MPI tests), then
# CTest will run 10 tests at the same time and will launch new tests as
# running tests finish.  One can also define tests which use more than one
# process or use more cores than the number of MPI processes.  When passing in
# ``NUM_MPI_PROCS <numMpiProcs>`` (see above), this TriBITS function will set the
# built-in CTest property ``PROCESSORS`` to ``<numMpiProcs>`` using::
#
#   set_tests_properties(<fullTestName> PROPERTIES PROCESSORS <numMpiProcs>)
#
# This tells CTest that the defined test uses ``<numMpiProcs>`` processes and
# CTest will use that information to not exceed the requested parallel level.
# For example, if several ``NUM_MPI_PROCS 3`` tests are defined and CTest is
# run with ``'ctest -j12'``, then CTest would schedule and run 4 of these
# tests at a time (to make use of 12 processors/cores on the machine),
# starting new tests as running tests finish, until all of the tests have been
# run.
#
# There are some situations where a test will use more processes/cores than
# specified by ``NUM_MPI_PROCS <numMpiProcs>`` such as when the underlying
# executable fires off more processes in parallel to do processing.  Also, an
# MPI program may use threading and therefore use overall more cores than the
# number of MPI processes. For these cases, it is critical to set
# ``NUM_TOTAL_CORES_USED <numTotalCoresUsed>`` to tell TriBITS and CTest how
# many cores will be used by a threaded test.  This is needed to exclude the
# test if there are too many processes/cores needed to run the test than are
# available.  Also, if the test is added, then this is needed to set the
# built-in CTest ``PROCESSORS`` property so CTest can avoid overloading the
# machine.  For example, for test where the MPI executable running on 4
# processes uses 10 threads per process, one would set::
#
#    NUM_MPI_PROCS  4  NUM_TOTAL_CORES_USED  40
#
# In this case, it sets the CTest ``PROCESSORS`` property as::
#
#   set_tests_properties(<fullTestName> PROPERTIES PROCESSORS <numTotalCoresUsed>)
#
# When the number of processes a test uses does not cleanly divide into the
# requested CTest parallel level, it is not clear how CTest schedules the
# tests (hard to find documentation on this but one could always inspect the
# CTest source code to find out for sure).  However, one boundary case that is
# well observed is that CTest will run all defined tests regardless of the
# size of the ``PROCESSORS`` property or the value of
# ``CTEST_PARALLEL_LEVEL``.  For example, if there are tests where
# ``PROCESSORS`` is set to 20 but ```ctest -j10'`` is used, then CTest will
# still run those tests (using 20 processes) but will not schedule any other
# tests while the parallel level is exceeded.  This can overload the machine
# obviously.  Therefore, always set ``MPI_EXEC_MAX_NUMPROCS`` to the maximum
# number of cores/processes that can be comfortably run on a given machine.
# Also note that MPI tests are very fragile to overloaded machines
#
# NOTE: **Never** manually override the ``PROCESSORS`` CTest property.
# Instead, always using ``NUM_TOTAL_CORES_USED <numTotalCoresUsed>`` to set
# this.  This is important because TriBITS needs to know how many
# processes/cores are required for test so that it can disable the test if the
# test requires more cores/processes than a given machine can handle or to
# exceed an imposed budget of the number processes to be used.  The latter is
# important when running multiple ``ctest -J<N>`` invocations on the same test
# machine.
#
# .. _Setting timeouts for tests (tribits_add_test()):
#
# **Setting timeouts for tests (tribits_add_test())**
#
# By default, all tests have a default timeout (1500 seconds for most
# projects, see `DART_TESTING_TIMEOUT`_).  That means that if they hang
# (e.g. as is common when deadlocking occurs with multi-process MPI-based
# tests and multi-threaded tests) then each test may hang for a long time,
# causing the overall test suite to take a long time to complete.  For many
# CMake projects, this default timeout is way too long.
#
# Timeouts for tests can be adjusted in a couple of ways.  First, a default
# timeout for all tests is enforced by CTest given the configured value of the
# variable `DART_TESTING_TIMEOUT`_ (typically set by the user but set by
# default to 1500 for most projects).  This is a global timeout that applies
# to all tests that don't otherwise have individual timeouts set using the
# ``TIMEOUT`` CTest property (see below).  The value of
# ``DART_TESTING_TIMEOUT`` in the CMake cache on input to CMake will get
# scaled by `${PROJECT_NAME}_SCALE_TEST_TIMEOUT`_ and the scaled value gets
# written into the file ``DartConfiguration.tcl`` as the field ``TimeOut``.
# The``ctest`` executable reads ``TimeOut`` from this file when it runs to
# determine the default global timeout.  The value of this default global
# ``TimeOut`` can be overridden using the ``ctest`` argument ``--timeout
# <seconds>`` (see `Overriding test timeouts`_).
#
# Alternatively, timeouts for individual tests can be set using the input
# argument ``TIMEOUT`` (see `Formal Arguments (tribits_add_test())`_ above).
# The timeout value passed in to this function is then scaled by
# ``${PROJECT_NAME}_SCALE_TEST_TIMEOUT`` and the scaled timeout is then set as
# the CTest test property ``TIMEOUT``.  One can observe the value of this
# property in the CMake-generated file ``CTestTestfile.cmake`` in the current
# build directory.  Individual test timeouts set this way are not impacted by
# the global default timeout ``DART_TESTING_TIMEOUT`` or the ``ctest``
# argument ``--timeout <seconds>``.
#
# In summary, CTest determines the timeout for any individual test as follows:
#
# * If the ``TIMEOUT`` CTest property in ``CTestTestfile.cmake`` for that test
#   is set, then that value is used.
#
# * Else if the ``ctest`` commandline option ``--timeout <seconds>`` is
#   provided, then that timeout is used.
#
# * Else if the property ``TimeOut`` in the base file
#   ``DartConfiguration.tcl`` is set to non-empty, then that timeout is used.
#
# * Else, no timeout is set and the test can run (and hang) forever until
#   manually killed by the user.
#
# .. _Debugging and Examining Test Generation (tribits_add_test()):
#
# **Debugging and Examining Test Generation (tribits_add_test())**
#
# In order to see what tests get added and if not then why, configure with
# ``${PROJECT_NAME}_TRACE_ADD_TEST=ON``.  That will print one line per show
# that the test got added and if not then why the test was not added (i.e. due
# to ``COMM``, ``NUM_MPI_PROCS``, ``CATEGORIES``, ``HOST``, ``XHOST``,
# ``HOSTTYPE``, or ``XHOSTTYPE``).
#
# Also, CMake writes a file ``CTestTestfile.cmake`` in the current binary
# directory which contains all of the added tests and test properties that are
# set.  This is the file that is read by ``ctest`` when it runs to determine
# what tests to run, determine pass/fail and adjust other behavior using test
# properties.  In this file, one can see the exact ``add_test()`` and
# ``set_tests_properties()`` commands.  The is the ultimate way to debug
# exactly what tests are getting added by this function (or if the test is
# even being added at all).
#
# .. _Disabling Tests Externally (tribits_add_test()):
#
# **Disabling Tests Externally (tribits_add_test())**
#
# The test can be disabled externally by setting the CMake cache variable
# ``<fullTestName>_DISABLE=TRUE``.  This allows tests to be disabled on a
# case-by-case basis by the user (for whatever reason).  Here,
# ``<fullTestName>`` must be the *exact* name that shows up in 'ctest -N' when
# running the test.  If multiple tests are added in this function through
# multiple argument sets to ``ARGS`` or through multiple
# ``POSTFIX_AND_ARGS_<IDX>`` arguments, then ``<fullTestName>_DISABLE=TRUE``
# must be set for each test individually.  When a test is disabled in this
# way, TriBITS will always print a warning to the ``cmake`` stdout at
# configure time warning that the test is being disabled.
#
# .. _Adding extra commandline arguments externally (tribits_add_test()):
#
# **Adding extra commandline arguments externally (tribits_add_test())**
#
# One can add additional command-line arguments for any ctest test added using
# this function.  In order to do so, set the CMake cache variable::
#
#   set(<fullTestName>_EXTRA_ARGS  "<earg0>;<earg1>;<earg2>;..."
#     CACHE  STRING  "Extra args")
#
# in a ``*.cmake`` configure options fragment file or::
#
#   -D <fullTestName>_EXTRA_ARGS="<earg0>;<earg1>;<earg2>;..."
#
# on the CMake command-line.
#
# These extra command-line arguments are added after any arguments passed in
# through ``ARGS "<oarg0> <oarg1> ..."`` or ``POSTFIX_AND_ARGS_<IDX> <oarg0>
# <oarg1> ...``.  This allows these extra arguments to override the earlier
# arguments.
#
# The primary motivating use case for ``<fullTestName>_EXTRA_ARGS`` is to
# allow one to alter how a test runs on a specific platform or build.  For
# example, this allows one to disable specific individual unit tests for a
# GTest executable such as with::
#
#   set(<fullTestName>_EXTRA_ARGS "--gtest_filter=-<unittest0>:<unittest1>:..."
#     CACHE  STRING  "Disable specific unit tests" )  
#
# For example, this would be an alternative to disabling an entire unit
# testing executable using ``-D<fullTestName>_DISABLE=ON`` as described above.
#
function(tribits_add_test EXE_NAME)

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("")
    message("TRIBITS_ADD_TEST: ${EXE_NAME} ${ARGN}")
  endif()

  global_set(TRIBITS_ADD_TEST_ADD_TEST_INPUT)
  global_set(TRIBITS_SET_TEST_PROPERTIES_INPUT)
  global_set(MESSAGE_WRAPPER_INPUT)

  #
  # A) Parse the input arguments
  #

  # Allow for a maximum of 20 (0 through 19) postfix and argument blocks
  set(MAX_NUM_POSTFIX_AND_ARGS_IDX 19)

  set(POSTFIX_AND_ARGS_LIST "")
  foreach( POSTFIX_AND_ARGS_IDX RANGE ${MAX_NUM_POSTFIX_AND_ARGS_IDX})
    list( APPEND POSTFIX_AND_ARGS_LIST POSTFIX_AND_ARGS_${POSTFIX_AND_ARGS_IDX} )
  endforeach()
  #print_var(POSTFIX_AND_ARGS_LIST)

  cmake_parse_arguments(
     PARSE_ARGV 1 # One named argument
     #prefix
     PARSE
     # options
     "NOEXEPREFIX;NOEXESUFFIX;STANDARD_PASS_OUTPUT;WILL_FAIL;ADD_DIR_TO_NAME;RUN_SERIAL"
     #one_value_keywords
     "DISABLED;LIST_SEPARATOR"
     #multi_value_keywords
"DIRECTORY;KEYWORDS;COMM;NUM_MPI_PROCS;NUM_TOTAL_CORES_USED;ARGS;${POSTFIX_AND_ARGS_LIST};NAME;NAME_POSTFIX;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;EXCLUDE_IF_NOT_TRUE;PASS_REGULAR_EXPRESSION;FAIL_REGULAR_EXPRESSION;TIMEOUT;ENVIRONMENT;ADDED_TESTS_NAMES_OUT"
     )

  tribits_check_for_unparsed_arguments()

  if (NOT "${PARSE_ARGS}" STREQUAL "")
    list(LENGTH PARSE_ARGS NUM_PARSE_ARGS)
  elseif (PARSE_POSTFIX_AND_ARGS_0)
    # We will use this list instead
  else()
    # Neither 'ARGS' nor 'POSTFIX_AND_ARGS' was selected so just assume one
    # empty arg
    set(PARSE_ARGS " ")
    set(NUM_PARSE_ARGS 1)
  endif()

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("")
    message("TRIBITS_ADD_TEST: EXE_NAME = ${EXE_NAME}")
  endif()

  if(PARSE_ADDED_TESTS_NAMES_OUT)
    set(${PARSE_ADDED_TESTS_NAMES_OUT} "" PARENT_SCOPE )
  endif()

  set(ADDED_TESTS_NAMES_OUT)

  #
  # Get test name
  #

  # If requested create a modifier for the name that will be inserted between
  # the package name and the given name or exe_name for the test
  set(directoryName "")
  if (PARSE_ADD_DIR_TO_NAME)
    tribits_create_name_from_current_source_directory(directoryName)
    set(directoryName "${directoryName}_")
  endif()

  #message("TRIBITS_ADD_TEST: ${EXE_NAME}: EXE_BINARY_NAME = ${EXE_BINARY_NAME}")

  if (PARSE_NAME)
    set(TEST_NAME "${directoryName}${PARSE_NAME}")
  elseif (PARSE_NAME_POSTFIX)
    set(TEST_NAME "${directoryName}${EXE_NAME}_${PARSE_NAME_POSTFIX}")
  else()
    set(TEST_NAME "${directoryName}${EXE_NAME}")
  endif()

  set(TEST_NAME "${PACKAGE_NAME}_${TEST_NAME}")

  #
  # B) Add or don't add tests based on a number of criteria
  #

  set(ADD_THE_TEST FALSE)
  tribits_add_test_process_enable_tests(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    return()
  endif()

  tribits_add_test_process_skip_ctest_add_test(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    return()
  endif()

  set(ADD_THE_TEST FALSE)
  tribits_add_test_process_categories(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    return()
  endif()

  set(ADD_THE_TEST FALSE)
  tribits_add_test_process_host_hosttype(ADD_THE_TEST)
  if (NOT ADD_THE_TEST)
    return()
  endif()

  #
  # C) Set the name and path of the binary that will be run
  #

  tribits_add_test_get_exe_binary_name(
    "${EXE_NAME}"
    ${PARSE_NOEXEPREFIX}
    ${PARSE_NOEXESUFFIX}
    ${PARSE_ADD_DIR_TO_NAME} EXE_BINARY_NAME
    )

  tribits_add_test_adjust_directory( ${EXE_BINARY_NAME} "${PARSE_DIRECTORY}"
    EXECUTABLE_PATH)

  #message("TRIBITS_ADD_TEST: ${EXE_NAME}: EXECUTABLE_PATH = ${EXECUTABLE_PATH}")

  #
  # D) Determine if we will add the serial or MPI tests based on input COMM
  # and TPL_ENABLE_MPI
  #

  tribits_process_comm_args(ADD_SERIAL_TEST  ADD_MPI_TEST  ${PARSE_COMM})

  #
  # E) Get the MPI options
  #

  tribits_add_test_get_num_procs_used("${PARSE_NUM_MPI_PROCS}"
    "NUM_MPI_PROCS"  NUM_PROCS_USED  NUM_PROCS_USED_NAME)
  if (NUM_PROCS_USED LESS 0)
    set(ADD_SERIAL_TEST FALSE)
    set(ADD_MPI_TEST FALSE)
  endif()

  if (TPL_ENABLE_MPI)
    set(MPI_NAME_POSTFIX "_MPI_${NUM_PROCS_USED}")
  else()
    set(MPI_NAME_POSTFIX "")
  endif()

  tribits_add_test_get_num_total_cores_used("${TEST_NAME}${MPI_NAME_POSTFIX}"
    "${PARSE_NUM_TOTAL_CORES_USED}"  "NUM_TOTAL_CORES_USED"
    "${NUM_PROCS_USED}"  "${NUM_PROCS_USED_NAME}"
    NUM_TOTAL_CORES_USED  SKIP_TEST)
  if (SKIP_TEST)
    set(ADD_SERIAL_TEST FALSE)
    set(ADD_MPI_TEST FALSE)
  endif()

  #
  # F) Add the tests
  #

  if (NOT ADD_SERIAL_TEST AND NOT ADD_MPI_TEST)
    return()
  endif()

  if (NOT "${PARSE_ARGS}" STREQUAL "")

    # F.1) Add tests with simple lists of arguments

    set(COUNTER 0)

    foreach(PARSE_ARG ${PARSE_ARGS})

      if(${NUM_PARSE_ARGS} EQUAL 1)
        set(TEST_NAME_INSTANCE "${TEST_NAME}${MPI_NAME_POSTFIX}")
      else()
        set(TEST_NAME_INSTANCE "${TEST_NAME}_${COUNTER}${MPI_NAME_POSTFIX}")
      endif()
      if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message(STATUS "TEST_NAME = ${TEST_NAME_INSTANCE}")
      endif()

      tribits_convert_cmnd_arg_string_to_add_test_arg_array(${PARSE_ARG} INARGS)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        print_var(INARGS)
      endif()

      tribits_set_run_serial(${TEST_NAME_INSTANCE} "${PARSE_RUN_SERIAL}"
        SET_RUN_SERIAL)

      tribits_set_disabled_and_msg(${TEST_NAME_INSTANCE} "${PARSE_DISABLED}"
        SET_DISABLED_AND_MSG)

      tribits_add_test_add_test_all( ${TEST_NAME_INSTANCE}
        "${EXECUTABLE_PATH}" "${PARSE_CATEGORIES}"  "${NUM_PROCS_USED}"
        "${NUM_TOTAL_CORES_USED}"
        "${SET_RUN_SERIAL}" "${SET_DISABLED_AND_MSG}" ADDED_TEST_NAME  ${INARGS}
        "${${TEST_NAME_INSTANCE}_EXTRA_ARGS}" )
      if(PARSE_ADDED_TESTS_NAMES_OUT AND ADDED_TEST_NAME)
        list(APPEND ADDED_TESTS_NAMES_OUT ${ADDED_TEST_NAME})
      endif()

      math(EXPR COUNTER ${COUNTER}+1 )

    endforeach()

  elseif (PARSE_POSTFIX_AND_ARGS_0)

    # F.2) Add tests with different postfixes for each set of arguments

    foreach( POSTFIX_AND_ARGS_IDX RANGE ${MAX_NUM_POSTFIX_AND_ARGS_IDX})

      if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
        print_var(PARSE_POSTFIX_AND_ARGS_${POSTFIX_AND_ARGS_IDX})
      endif()

      if (NOT PARSE_POSTFIX_AND_ARGS_${POSTFIX_AND_ARGS_IDX})
        break()
      endif()

      set( POSTFIX_AND_ARGS "${PARSE_POSTFIX_AND_ARGS_${POSTFIX_AND_ARGS_IDX}}" )

      list( GET  POSTFIX_AND_ARGS  0  POSTFIX )
      set( INARGS  "${POSTFIX_AND_ARGS}" ) # Initially contains postfix as ele 0
      list( REMOVE_AT  INARGS  0 ) # Strip off the postfix name

      if (NOT "${POSTFIX}" STREQUAL "")
        set(TEST_NAME_INSTANCE "${TEST_NAME}_${POSTFIX}${MPI_NAME_POSTFIX}")
      else()
        set(TEST_NAME_INSTANCE "${TEST_NAME}${MPI_NAME_POSTFIX}")
      endif()

      tribits_set_run_serial(${TEST_NAME_INSTANCE} "${PARSE_RUN_SERIAL}"
        SET_RUN_SERIAL)

      tribits_set_disabled_and_msg(${TEST_NAME_INSTANCE} "${PARSE_DISABLED}"
        SET_DISABLED_AND_MSG)

      tribits_add_test_add_test_all( ${TEST_NAME_INSTANCE}
        "${EXECUTABLE_PATH}" "${PARSE_CATEGORIES}" "${NUM_PROCS_USED}" 
        "${NUM_TOTAL_CORES_USED}"
        "${SET_RUN_SERIAL}" "${SET_DISABLED_AND_MSG}" ADDED_TEST_NAME  ${INARGS}
        "${${TEST_NAME_INSTANCE}_EXTRA_ARGS}"
        )
      if(PARSE_ADDED_TESTS_NAMES_OUT AND ADDED_TEST_NAME)
        list(APPEND ADDED_TESTS_NAMES_OUT ${ADDED_TEST_NAME})
      endif()

    endforeach()

  endif()

  if(PARSE_ADDED_TESTS_NAMES_OUT)
    set(${PARSE_ADDED_TESTS_NAMES_OUT} "${ADDED_TESTS_NAMES_OUT}"
      PARENT_SCOPE )
  endif()

endfunction()
