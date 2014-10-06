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

INCLUDE(TribitsAddTestHelpers)


#
# @FUNCTION: TRIBITS_ADD_TEST()
#
# Add a test or a set of tests for a single executable or command using CTest
# ``ADD_TEST()``.
#
# Usage::
#
#   TRIBITS_ADD_TEST(
#     <exeRootName>  [NOEXEPREFIX]  [NOEXESUFFIX]
#     [NAME <testName> | NAME_POSTFIX <testNamePostfix>]
#     [DIRECTORY <directory>]
#     [ADD_DIR_TO_NAME]
#     [ARGS "<arg0> <arg1> ..." "<arg2> <arg3> ..." ...
#       | POSTFIX_AND_ARGS_0 <postfix0> <arg0> <arg1> ...
#         POSTFIX_AND_ARGS_1 ... ]
#     [COMM [serial] [mpi]]
#     [NUM_MPI_PROCS <numProcs>]
#     [CATEGORIES <category0>  <category1> ...]
#     [HOST <host0> <host1> ...]
#     [XHOST <host0> <host1> ...]
#     [HOSTTYPE <hosttype0> <hosttype1> ...]
#     [XHOSTTYPE <hosttype0> <hosttype1> ...]
#     [STANDARD_PASS_OUTPUT
#       | PASS_REGULAR_EXPRESSION "<regex0>;<regex1>;..."]
#     [FAIL_REGULAR_EXPRESSION "<regex0>;<regex1>;..."]
#     [WILL_FAIL]
#     [ENVIRONMENT <var0>=<value0> <var1>=<value1> ...]
#     [TIMEOUT <maxSeconds>]
#     )
#
# *Sections:*
#
# * `Formal Arguments (TRIBITS_ADD_TEST())`_
# * `Determining the Executable or Command to Run (TRIBITS_ADD_TEST())`_
# * `Determining the Full Test Name (TRIBITS_ADD_TEST())`_
# * `Adding Multiple Tests  (TRIBITS_ADD_TEST())`_
# * `Determining Pass/Fail (TRIBITS_ADD_TEST())`_
# * `Setting additional test properties (TRIBITS_ADD_TEST())`_
# * `Running multiple tests at the same time (TRIBITS_ADD_TEST())`_
# * `Debugging and Examining Test Generation (TRIBITS_ADD_TEST())`_
# * `Disabling Tests Externally (TRIBITS_ADD_TEST())`_
#
# .. _Formal Arguments (TRIBITS_ADD_TEST()):
#
# **Formal Arguments (TRIBITS_ADD_TEST())**
#
#   ``<exeRootName>``
#
#     The name of the executable or path to the executable to run for the test
#     (see `Determining the Executable or Command to Run
#     (TRIBITS_ADD_TEST())`_).  This name is also the default root name for
#     the test (see `Determining the Full Test Name (TRIBITS_ADD_TEST())`_).
#
#   ``NOEXEPREFIX``
#
#    If specified, then the prefix ``${PACKAGE_NAME}_`` is assumed **not** to
#    be prepended to ``<exeRootName>`` (see `Determining the Executable or
#    Command to Run (TRIBITS_ADD_TEST())`_).
#
#   ``NOEXESUFFIX``
#
#      If specified, then the postfix
#      ``${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}`` is assumed **not** to be
#      post-pended to ``<exeRootName>`` (see `Determining the Executable or
#      Command to Run (TRIBITS_ADD_TEST())`_).
#
#   ``NAME <testRootName>``
#
#     If specified, gives the root name of the test.  If not specified, then
#     ``<testRootName>`` is taken to be ``<exeRootName>``.  The actual test
#     name passed to ``ADD_TEST()`` will always be prefixed as
#     ``${PACKAGE_NAME}_<testRootName>``.  The main purpose of this argument
#     is to allow multiple tests to be defined for the same executable.  CTest
#     requires all test names to be globally unique in a single project.  See
#     `Determining the Full Test Name (TRIBITS_ADD_TEST())`_.
#
#   ``NAME_POSTFIX <testNamePostfix>``
#
#     If specified, gives a postfix that will be added to the standard test
#     name based on ``<exeRootName>`` (appended as ``_<NAME_POSTFIX>``).  If
#     the ``NAME <testRootName>`` argument is given, this argument is ignored.
#     See `Determining the Full Test Name (TRIBITS_ADD_TEST())`_.
#
#   ``DIRECTORY <dir>``
#
#     If specified, then the executable is assumed to be in the directory
#     given by ``<dir>``.  The directory ``<dir>`` can either be a relative or
#     absolute path.  If not specified, the executable is assumed to be in the
#     current binary directory ``${CMAKE_CURRENT_BINARY_DIR}``.  See
#     `Determining the Executable or Command to Run (TRIBITS_ADD_TEST())`_.
#
#   ``ADD_DIR_TO_NAME``
#
#     If specified, then the directory name that this test resides in will be
#     added into the name of the test after the package name is added and
#     before the root test name (see `Determining the Full Test Name
#     (TRIBITS_ADD_TEST())`_).  The directory name will have the package's
#     base directory stripped off so only the unique part of the test
#     directory will be used.  All directory separators ``"/"`` will be
#     changed into underscores ``"_"``.
#
#   ``RUN_SERIAL``
#
#     If specified then no other tests will be allowed to run while this test
#     is running. This is useful for devices (like CUDA GPUs) that require
#     exclusive access for processes/threads.  This just sets the CTest test
#     property ``RUN_SERIAL`` using the built-in CMake function
#     ``SET_TESTS_PROPERTIES()``.
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
#     Multiple Tests (TRIBITS_ADD_TEST())`_ for more details and exmaples.
#
#   ``POSTFIX_AND_ARGS_<IDX> <postfix> <arg0> <arg1> ...``
#
#     If specified, gives a sequence of sets of test postfix names and
#     arguments lists for different tests (up to ``POSTFIX_AND_ARGS_19``).
#     For example, a set of three different tests with argument lists can be
#     specified as::
#
#       POSTIFX_AND_ARGS_0 postfix0 --arg1 --arg2="dummy"
#       POSTIFX_AND_ARGS_1 postfix1  --arg2="fly"
#       POSTIFX_AND_ARGS_2 postfix2  --arg2="bags"
#
#     This will create three different test cases with the postfix names
#     ``postfix0``, ``postfix1``, and ``postfix2``.  The indexes must be
#     consecutive starting a ``0`` and going up to (currently) ``19``.  The
#     main advantages of using these arguments instead of just ``ARGS`` are
#     that one can give a meaningful name to each test case and one can
#     specify multiple arguments without having to quote them and one can
#     allow long argument lists to span multiple lines.  See `Adding Multiple
#     Tests (TRIBITS_ADD_TEST())`_ for more details and exmaples.
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
#   ``NUM_MPI_PROCS <numProcs>``
#
#     If specified, gives the number of MPI processes used to run the test
#     with the MPI exec program ``${MPI_EXEC}``.  If ``<numProcs>`` is greater
#     than ``${MPI_EXEC_MAX_NUMPROCS}`` then the test will be excluded.  If
#     not specified, then the default number of processes for an MPI build
#     (i.e. ``TPL_ENABLE_MPI=ON``) will be ``${MPI_EXEC_DEFAULT_NUMPROCS}``.
#     For serial builds (i.e. ``TPL_ENABLE_MPI=OFF``), this argument is
#     ignored.  This will also be set as the built-in test property
#     ``PROCESSORS`` to tell CTest how many processes this test will use (see
#     `Running multiple tests at the same time (TRIBITS_ADD_TEST())`_).
#
#   ``CATEGORIES <category0> <category1> ...``
#
#     If specified, gives the specific categories of the test.  Valid test
#     categories include ``BASIC``, ``CONTINUOUS``, ``NIGHTLY``, ``WEEKLY``
#     and ``PERFORMANCE``.  If not specified, the default category is
#     ``BASIC``.  When the test category does not match
#     ``${PROJECT_NAME}_TEST_CATEGORIES``, then the test is **not** added.
#     When ``CATEGORIES`` contains ``BASIC`` it will match
#     ``${PROJECT_NAME}_TEST_CATEGORIES`` equal to ``CONTINUOUS``,
#     ``NIGHTLY``, and ``WEEKLY``.  When ``CATEGORIES`` contains
#     ``CONTINUOUS`` it will match ``${PROJECT_NAME}_TEST_CATEGORIES`` equal
#     to ``CONTINUOUS``, ``NIGHTLY``, and ``WEEKLY``.  When ``CATEGORIES``
#     contains ``NIGHTLY`` it will match ``${PROJECT_NAME}_TEST_CATEGORIES``
#     equal to ``NIGHTLY`` and ``WEEKLY``.  When ``CATEGORIES`` contains
#     ``PERFORMANCE`` it will match
#     ``${PROJECT_NAME}_TEST_CATEGORIES=PERFORMANCE`` only.
#
#   ``HOST <host0> <host1> ...``
#
#     If specified, gives a list of hostnames where the test will be included.
#     The current hostname is determined by the built-in CMake command
#     ``SITE_NAME(${PROJECT_NAME}_HOSTNAME)``.  On Linux/Unix systems, this is
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
#     added.  Typical host system type names include ``Linux``, ``Darwain``,
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
#   ``STANDARD_PASS_OUTPUT``
#
#     If specified, then the standard test output string ``End Result: TEST
#     PASSED`` is grepped in the test stdout for to determine success.  This
#     is needed for MPI tests on some platforms since the return value from
#     MPI executables is unreliable.  This is set using the built-in CTest
#     property ``PASS_REGULAR_EXPRESSION``.
#
#   ``PASS_REGULAR_EXPRESSION "<regex0>;<regex1>;..."``
#
#     If specified, then the test will be assumed to pass only if one of the
#     regular expressions ``<regex0>``, ``<regex1>`` etc. match the output
#     send to stdout.  Otherwise, the test will fail.  This is set using the
#     built-in CTest property ``PASS_REGULAR_EXPRESSION``.  Consult standard
#     CMake documentation for full behavior.
#
#   ``FAIL_REGULAR_EXPRESSION "<regex0>;<regex1>;..."``
#
#     If specified, then a test will be assumed to fail if one of the regular
#     expressions ``<regex0>``, ``<regex1>`` etc. match the output send to
#     stdout.  Otherwise, the test will pass.  This is set using the built-in
#     CTest property ``FAIL_REGULAR_EXPRESSION``.  Consult standard CMake
#     documentation for full behavior.
#
#   ``WILL_FAIL``
#
#     If passed in, then the pass/fail criteria will be inverted.  This is set
#     using the built-in CTest property ``WILL_FAIL``.  Consult standard CMake
#     documentation for full behavior.
#
#   ``ENVIRONMENT <var0>=<value0> <var1>=<value1> ...``
#
#     If passed in, the listed environment variables will be set before
#     calling the test.  This is set using the built-in CTest property
#     ``ENVIRONMENT``.
#
#   ``TIMEOUT <maxSeconds>``
#
#     If passed in, gives maximum number of seconds the test will be allowed
#     to run before being timed-out.  This sets the CTest property
#     ``TIMEOUT``.  The value ``<maxSeconds>`` will be scaled by the value of
#     `${PROJECT_NAME}_SCALE_TEST_TIMEOUT`_.
#
#     **WARNING:** Rather than just increasing the timeout for an expensive
#     test, please try to either make the test run faster or relegate the test
#     to being run less often (i.e. set ``CATEGORIES NIGHTLY`` or even
#     ``WEEKLY`` for extremely expensive tests).  Expensive tests are one of
#     the worse forms of technical debt that a project can have!
#
# In the end, this function just calls the built-in CMake commands
# ``ADD_TEST(${TEST_NAME} ...)`` and ``SET_TESTS_PROPERTIES(${TEST_NAME}
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
# function `TRIBITS_ADD_ADVANCED_TEST()`_ should be used instead.
#
# In the following subsections, more details on how tests are defined and run
# is given.
#
# .. _Determining the Executable or Command to Run (TRIBITS_ADD_TEST()):
#
# **Determining the Executable or Command to Run (TRIBITS_ADD_TEST())**
#
# This function is primarily designed to make it easy to run tests for
# executables built using the function `TRIBITS_ADD_EXECUTABLE()`_.  To set up
# tests to run arbitrary executables, see below.
#
# By default, the executable to run is determined by first getting the
# executable name which by default is assumed to be::
#
#  <fullExeName> =
#    ${PACKAGE_NAME}_<exeRootName>${${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX}
#
# which is (by no coincidence) identical to how it is selected in
# `TRIBITS_ADD_EXECUTABLE()`_.  This name can be altered by passing in
# ``NOEXEPREFIX``, ``NOEXESUFFIX``, and ``ADD_DIR_TO_NAME`` as described in
# `Executable and Target Name (TRIBITS_ADD_EXECUTABLE())`_.
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
# Whatever executable path is specified using this logic, if the executable is
# not found, then when ``ctest`` goes to run the test, it will mark it as
# ``NOT RUN``.
#
# .. _Determining the Full Test Name (TRIBITS_ADD_TEST()):
#
# **Determining the Full Test Name (TRIBITS_ADD_TEST())**
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
# ``ARGS``, then this is the test name actually passed in to ``ADD_TEST()``.
# If multiple tests are defined, then this name becomes the base test name for
# each of the tests (see `Adding Multiple Tests (TRIBITS_ADD_TEST())`_).
#
# Finally, for any test that gets defined, if MPI is enabled
# (i.e. ``TPL_ENABLE_MPI=ON``), then the terminal suffix
# ``_MPI_${NUM_MPI_PROCS}`` will be added to the end of the test name (even
# for multiple tests).  No such prefix is added for the serial case
# (i.e. ``TPL_ENABLE_MPI=OFF``).
#
# .. _Adding Multiple Tests  (TRIBITS_ADD_TEST()):
#
# **Adding Multiple Tests  (TRIBITS_ADD_TEST())**
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
# quoted set of arugments such as with::
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
# arugments ``<arg0>``, ``<arg1>``, ... do not need to be quoted and can
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
# .. _Determining Pass/Fail (TRIBITS_ADD_TEST()):
#
# **Determining Pass/Fail (TRIBITS_ADD_TEST())**
#
# The only means to determine pass/fail is to use the built-in CTest
# properties ``PASS_REGULAR_EXPRESSION`` and ``FAIL_REGULAR_EXPRESSION`` which
# can only grep the test's STDOUT/STDERR or to check for a 0 return value (or
# invert these using ``WILL_FAIL``).  For simple tests, that is enough.
# However, for more complex executables, one may need to examine one or more
# output files to determine pass/fail.  Raw CMake/CTest cannot do this.  In
# this case, one should use `TRIBITS_ADD_ADVANCED_TEST()`_ instead to add the
# test.
#
# .. _Setting additional test properties (TRIBITS_ADD_TEST()):
#
# **Setting additional test properties (TRIBITS_ADD_TEST())**
#
# After this function returns, any tests that get added using ``ADD_TEST()``
# can have additional properties set and changed using
# ``SET_TESTS_PROPERTIES()``.  Therefore, any tests properties that are not
# directly supported and passed through this wrapper function can be set in
# the outer ``CMakeLists.txt`` file after the call to ``TRIBITS_ADD_TEST()``.
#
# ToDo: Describe how to use new variable ``ADDED_TESTS_OUT`` to get the list
# of tests actually added (if they are added) in order to make it easy to set
# additional test properties.
#
# .. _Running multiple tests at the same time (TRIBITS_ADD_TEST()):
#
# **Running multiple tests at the same time (TRIBITS_ADD_TEST())**
#
# By default, CTest will run many tests defined with ``ADD_TEST()`` at same
# time as it can according to its parallel level (e.g. ``'test -j<N>'`` or the
# CTest property ``CTEST_PARALLEL_LEVEL``).  For example, when raw ``'ctest
# -j10'`` is run, CTest will run multiple tests at the same time to try to
# make usage of 10 processes.  If all of the defined tests only used one
# process (which is assumed by default except for MPI tests), then CTest will
# run 10 tests at the same time and will launch new tests as running tests
# finish.  One can also define tests using ``ADD_TEST()`` that use more than
# one process such as for MPI tests.  When passing in ``NUM_MPI_PROCS
# <numProcs>`` (see above), this TriBITS function will set the built-in CTest
# property ``PROCESSORS`` to ``<numProcs>`` using::
#
#   SET_TESTS_PROPERTIES(<fullTestName> PROPERTIES PROCESSORS <numProcs>)
#
# This tells CTest that the defined test uses ``<numProcs>`` processes and
# CTest will use that information to not exceed the requested parallel level.
# For example, if several ``NUM_MPI_PROCS 3`` tests are defined and CTest is
# run with ``'ctest -j12'``, then CTest would schedule and run 4 of these
# tests at a time (to make use of 12 processes), starting new ones as running
# tests finish, until all of the tests have been run.
#
# When the number of processes a test uses does not cleanly divide into the
# requested CTest parallel level, it is not clear how CTest schedules the
# tests (hard to find documentation on this but one could always inspect the
# CTest source code to find out for sure).  However, one boundary case that is
# well observed is that CTest will run all defined tests regardless of the
# size of the ``PROCESSORS`` property or the value of
# ``CTEST_PARALLEL_LEVEL``.  For example, if there are tests where
# ``PROCESSORS`` is set to 20 but ```ctest -j10'`` is run, then CTest will
# still run those tests (using 20 processes) one at a time but will not
# schedule any other tests while the parallel level is exceeded.
#
# For single-thread MPI tests, the behavior built into TriBITS does exactly
# the right thing.  Defining the test with ``NUM_MPI_PROCS <numProcs>`` will
# call ``${MPI_EXEC}`` with ``<numProcs>`` and it will set the CTest property
# ``PROCESSORS`` to ``<numProcs>``.  However, if the MPI processes use more
# than one thread, then CTest could easily oversubscribe the machine.  For
# example, consider the case where one is on a machine that only has 16 cores
# and one defines MPI tests with ``NUM_MPI_PROCS 4`` but each MPI process
# launches 6 threads.  In this case, running these tests with ``'ctest -j8'``,
# CTest would schedule 2 of these 4-process tests to run at a time but would
# in actuality be using ``2*4*6 = 48`` cores and would overload 32 core
# machine.  The other case that is not automatically handled by TriBITS is
# when a test script (not MPI) launches multiple processes simultaneously
# internally.
#
# Therefore, in cases where the executable or script uses multiple processes,
# then one must manually override the ``PROCESSORS`` property.  To do, this
# after the ``TRIBITS_ADD_TEST()`` (or `TRIBITS_ADD_ADVANCED_TEST()`_)
# function returns, one can reset the ``PROCESSORS`` property` with::
#
#   SET_TESTS_PROPERTIES(<fullTestName> PROPERTIES PROCESSORS <fullNumProces>)
#
# For example, if one runs an MPI program that uses 4 processes and 6 threads
# per process, one would call::
#
#   TRIBITS_ADD_TEST(myProg ... NUM_MPI_PROCS 4 ...)
#   SET_TESTS_PROPERTIES(${PACKAGE_NAME}_myProg PROPERTIES PROCESSORS 12)
#
# ToDo: Update above example to use loop over ``ADDED_TESTS_OUT``.
#
# .. _Debugging and Examining Test Generation (TRIBITS_ADD_TEST()):
#
# **Debugging and Examining Test Generation (TRIBITS_ADD_TEST())**
#
# In order to see what tests are getting added and to debug some issues in
# test creation, one can set the cache variable
# ``${PROJECT_NAME}_VERBOSE_CONFIGURE=ON``.  This will result in the printout
# of some information about the test getting added or not.
#
# Also, CMake writes a file ``CTestTestfile.cmake`` in the current binary
# directory which contains all of the added tests and test properties that are
# set.  This is the file that is read by ``ctest`` when it runs to determine
# what tests to run, determine pass/fail and adjust other behavior using test
# properties.  In this file, one can see the exact ``ADD_TEST()`` and
# ``SET_TESTS_PROPERTIES()`` commands.  The is the ultimate way to debug
# exactly what tests are getting added by this function (or if the test is
# even being added at all).
#
# .. _Disabling Tests Externally (TRIBITS_ADD_TEST()):
#
# **Disabling Tests Externally (TRIBITS_ADD_TEST())**
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
FUNCTION(TRIBITS_ADD_TEST EXE_NAME)

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRIBITS_ADD_TEST: ${EXE_NAME} ${ARGN}")
  ENDIF()

  GLOBAL_SET(TRIBITS_ADD_TEST_ADD_TEST_INPUT "")
  GLOBAL_SET(TRIBITS_SET_TEST_PROPERTIES_INPUT)

  #
  # A) Parse the input arguments
  #

  # Allow for a maximum of 20 (0 through 19) postfix and argument blocks
  SET(MAX_NUM_POSTFIX_AND_ARGS_IDX 19)

  SET(POSTFIX_AND_ARGS_LIST "")
  FOREACH( POSTFIX_AND_ARGS_IDX RANGE ${MAX_NUM_POSTFIX_AND_ARGS_IDX})
    LIST( APPEND POSTFIX_AND_ARGS_LIST POSTFIX_AND_ARGS_${POSTFIX_AND_ARGS_IDX} )
  ENDFOREACH()
  #PRINT_VAR(POSTFIX_AND_ARGS_LIST)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "DIRECTORY;KEYWORDS;COMM;NUM_MPI_PROCS;ARGS;${POSTFIX_AND_ARGS_LIST};NAME;NAME_POSTFIX;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;PASS_REGULAR_EXPRESSION;FAIL_REGULAR_EXPRESSION;TIMEOUT;ENVIRONMENT;ADDED_TESTS_NAMES_OUT"
     #options
     "NOEXEPREFIX;NOEXESUFFIX;STANDARD_PASS_OUTPUT;WILL_FAIL;ADD_DIR_TO_NAME;RUN_SERIAL"
     ${ARGN}
     )

  IF (PARSE_ARGS)
    LIST(LENGTH PARSE_ARGS NUM_PARSE_ARGS)
  ELSEIF (PARSE_POSTFIX_AND_ARGS_0)
    # We will use this list instead
  ELSE()
    # Niether 'ARGS' nor 'POSTFIX_AND_ARGS' was selected so just assume one
    # empty arg
    SET(PARSE_ARGS " ")
    SET(NUM_PARSE_ARGS 1)
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRIBITS_ADD_TEST: EXE_NAME = ${EXE_NAME}")
  ENDIF()

  IF(PARSE_ADDED_TESTS_NAMES_OUT)
    SET(${PARSE_ADDED_TESTS_NAMES_OUT} PARENT_SCOPE)
  ENDIF()

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

  #
  # C) Set the name and path of the binary that will be run
  #

  TRIBITS_ADD_TEST_GET_EXE_BINARY_NAME(
    "${EXE_NAME}"
    ${PARSE_NOEXEPREFIX}
    ${PARSE_NOEXESUFFIX}
    ${PARSE_ADD_DIR_TO_NAME} EXE_BINARY_NAME
    )

  # If requested create a modifier for the name that will be inserted between
  # the package name and the given name or exe_name for the test
  SET(DIRECTORY_NAME "")
  IF(PARSE_ADD_DIR_TO_NAME)
    TRIBITS_CREATE_NAME_FROM_CURRENT_SOURCE_DIRECTORY(DIRECTORY_NAME)
    SET(DIRECTORY_NAME "${DIRECTORY_NAME}_")
  ENDIF()

  #MESSAGE("TRIBITS_ADD_TEST: ${EXE_NAME}: EXE_BINARY_NAME = ${EXE_BINARY_NAME}")

  IF (PARSE_NAME)
    SET(TEST_NAME "${DIRECTORY_NAME}${PARSE_NAME}")
  ELSEIF (PARSE_NAME_POSTFIX)
    SET(TEST_NAME "${DIRECTORY_NAME}${EXE_NAME}_${PARSE_NAME_POSTFIX}")
  ELSE()
    SET(TEST_NAME "${DIRECTORY_NAME}${EXE_NAME}")
  ENDIF()

  TRIBITS_ADD_TEST_ADJUST_DIRECTORY( ${EXE_BINARY_NAME} "${PARSE_DIRECTORY}"
    EXECUTABLE_PATH)

  #MESSAGE("TRIBITS_ADD_TEST: ${EXE_NAME}: EXECUTABLE_PATH = ${EXECUTABLE_PATH}")

  #
  # D) Determine if we will add the serial or MPI tests based on input COMM
  # and TPL_ENABLE_MPI
  #

  TRIBITS_PROCESS_COMM_ARGS(ADD_SERIAL_TEST  ADD_MPI_TEST  ${PARSE_COMM})

  #
  # E) Get the MPI options
  #

  TRIBITS_ADD_TEST_GET_NUM_PROCS_USED("${PARSE_NUM_MPI_PROCS}" NUM_PROCS_USED)
  IF (NUM_PROCS_USED LESS 0)
    SET(ADD_MPI_TEST FALSE)
  ENDIF()

  #
  # F) Add the tests
  #

  IF (NOT ADD_SERIAL_TEST AND NOT ADD_MPI_TEST)
    RETURN()
  ENDIF()

  IF (TPL_ENABLE_MPI)
    SET(MPI_NAME_POSTFIX "_MPI_${NUM_PROCS_USED}")
  ELSE()
    SET(MPI_NAME_POSTFIX "")
  ENDIF()

  IF (PARSE_ARGS)

    # F.1) Add tests with simple lists of arguments

    SET(COUNTER 0)

    FOREACH(PARSE_ARG ${PARSE_ARGS})

      IF(${NUM_PARSE_ARGS} EQUAL 1)
        SET(TEST_NAME_INSTANCE "${TEST_NAME}${MPI_NAME_POSTFIX}")
      ELSE()
        SET(TEST_NAME_INSTANCE "${TEST_NAME}_${COUNTER}${MPI_NAME_POSTFIX}")
      ENDIF()
      IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "TEST_NAME = ${TEST_NAME_INSTANCE}")
      ENDIF()

      TRIBITS_CONVERT_CMND_ARG_STRING_TO_ADD_TEST_ARG_ARRAY(${PARSE_ARG} INARGS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(INARGS)
      ENDIF()

      TRIBITS_ADD_TEST_ADD_TEST_ALL( ${TEST_NAME_INSTANCE}
        "${EXECUTABLE_PATH}"  "${NUM_PROCS_USED}"
        ${PARSE_RUN_SERIAL} ${INARGS} )
      IF(PARSE_ADDED_TESTS_NAMES_OUT)
        LIST(APPEND ${PARSE_ADDED_TESTS_NAMES_OUT}
          ${TEST_NAME_INSTANCE}
          PARENT_SCOPE
          )
      ENDIF()

      MATH(EXPR COUNTER ${COUNTER}+1 )

    ENDFOREACH()

  ELSEIF (PARSE_POSTFIX_AND_ARGS_0)

    # F.2) Add tests with different postfixes for each set of arguments

    FOREACH( POSTFIX_AND_ARGS_IDX RANGE ${MAX_NUM_POSTFIX_AND_ARGS_IDX})

      IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(PARSE_POSTFIX_AND_ARGS_${POSTFIX_AND_ARGS_IDX})
      ENDIF()

      IF (NOT PARSE_POSTFIX_AND_ARGS_${POSTFIX_AND_ARGS_IDX})
        BREAK()
      ENDIF()

      SET( POSTFIX_AND_ARGS ${PARSE_POSTFIX_AND_ARGS_${POSTFIX_AND_ARGS_IDX}} )

      LIST( GET  POSTFIX_AND_ARGS  0  POSTFIX )
      SET( INARGS  ${POSTFIX_AND_ARGS} ) # Initially contains postfix as ele 0
      LIST( REMOVE_AT  INARGS  0 ) # Strip off the postfix name

      SET(TEST_NAME_INSTANCE "${TEST_NAME}_${POSTFIX}${MPI_NAME_POSTFIX}")

      TRIBITS_ADD_TEST_ADD_TEST_ALL( ${TEST_NAME_INSTANCE}
        "${EXECUTABLE_PATH}"  "${NUM_PROCS_USED}" ${PARSE_CREATE_WORKING_DIR}
        ${PARSE_RUN_SERIAL} ${INARGS} )
      IF(PARSE_ADDED_TESTS_NAMES_OUT)
        LIST(APPEND ${PARSE_ADDED_TESTS_NAMES_OUT}
          ${TEST_NAME_INSTANCE}
          PARENT_SCOPE
          )
      ENDIF()

    ENDFOREACH()

  ENDIF()

ENDFUNCTION()
