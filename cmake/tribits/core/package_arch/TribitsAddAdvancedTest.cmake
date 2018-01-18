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
INCLUDE(TribitsConstants)

INCLUDE(AppendStringVar)
INCLUDE(PrintVar)


#
# @FUNCTION: TRIBITS_ADD_ADVANCED_TEST()
#
# Function that creates an advanced test defined by stringing together one or
# more executable and/or command invocations that is run as a ``cmake -P``
# script with very flexible pass/fail criteria.
#
# Usage::
#
#   TRIBITS_ADD_ADVANCED_TEST(
#     <testNameBase>
#     TEST_0 (EXEC <execTarget0> | CMND <cmndExec0>) ...
#     [TEST_1 (EXEC <execTarget1> | CMND <cmndExec1>) ...]
#     ...
#     [TEST_N (EXEC <execTargetN> | CMND <cmndExecN>) ...]
#     [OVERALL_WORKING_DIRECTORY (<overallWorkingDir> | TEST_NAME)]
#     [SKIP_CLEAN_OVERALL_WORKING_DIRECTORY]
#     [FAIL_FAST]
#     [RUN_SERIAL]
#     [KEYWORDS <keyword1> <keyword2> ...]
#     [COMM [serial] [mpi]]
#     [OVERALL_NUM_MPI_PROCS <overallNumProcs>]
#     [OVERALL_NUM_TOTAL_CORES_USED <overallNumTotalCoresUsed>]
#     [CATEGORIES <category0> <category1> ...]
#     [HOST <host0> <host1> ...]
#     [XHOST <host0> <host1> ...]
#     [HOSTTYPE <hosttype0> <hosttype1> ...]
#     [XHOSTTYPE <hosttype0> <hosttype1> ...]
#     [EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...]
#     [FINAL_PASS_REGULAR_EXPRESSION <regex> |
#       FINAL_FAIL_REGULAR_EXPRESSION <regex>]
#     [ENVIRONMENT <var1>=<value1> <var2>=<value2> ...]
#     [TIMEOUT <maxSeconds>]
#     [ADDED_TEST_NAME_OUT <testName>]
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
# For more details on these arguments, see `TEST_<idx> EXEC/CMND Test Blocks
# and Arguments (TRIBITS_ADD_ADVANCED_TEST())`_.
#
# The most common type of an atomic test block ``TEST_<idx>`` runs a command
# as either a package-built executable or just any command.  An atomic test
# command block ``TEST_<idx>`` (i.e. ``TEST_0``, ``TEST_1``, ...) takes the
# form::
#
#   TEST_<idx>
#      (EXEC <exeRootName> [NOEXEPREFIX] [NOEXESUFFIX] [ADD_DIR_TO_NAME]
#             [DIRECTORY <dir>]
#         | CMND <cmndExec>)
#      [ARGS <arg1> <arg2> ... <argn>]
#      [MESSAGE "<message>"]
#      [WORKING_DIRECTORY <workingDir>]
#      [SKIP_CLEAN_WORKING_DIRECTORY]
#      [NUM_MPI_PROCS <numProcs>]
#      [NUM_TOTAL_CORES_USED <numTotalCoresUsed>]
#      [OUTPUT_FILE <outputFile>]
#      [NO_ECHO_OUTPUT]]
#      [PASS_ANY
#        | PASS_REGULAR_EXPRESSION "<regex>"
#        | PASS_REGULAR_EXPRESSION_ALL "<regex1>" "<regex2>" ... "<regexn>"
#        | STANDARD_PASS_OUTPUT ]
#      [FAIL_REGULAR_EXPRESSION "<regex>"]
#      [ALWAYS_FAIL_ON_NONZERO_RETURN | ALWAYS_FAIL_ON_ZERO_RETURN]
#      [WILL_FAIL]
#
# For more information on these arguments, see `TEST_<idx> EXEC/CMND Test
# Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST())`_.
#
# The other type of ``TEST_<idx>`` block supported is for copying files and
# takes the form::
#
#   TEST_<idx>
#     COPY_FILES_TO_TEST_DIR <file0> <file1> ... <filen>
#     [SOURCE_DIR <srcDir>]
#     [DEST_DIR <destDir>]
#
# This makes it easy to copy files from the source tree (or other location) to
# inside of the test directory (usually created with ``OVERALL_WORKING_DIR
# TEST_NAME``) so that tests can run in their own private working directory
# (and so these files get deleted and recopied each time the test runs).  This
# approach has several advantages:
#
# * One can modify the input files and then just run the test with ``ctest``
#   in an iterative manner (and not have to configure again when using
#   ``CONFIGURE_FILE( ... COPYONLY)`` or build again when using
#   ``TRIBITS_COPY_FILES_TO_BINARY_DIR()`` in order to copy files).
#
# * When using ``OVERALL_WORKING_DIR TEST_NAME``, the test directory gets blow
#   away every time before it runs and therefore any old files are deleted
#   before the test gets run again (which avoids the problem of having a test
#   pass looking for the old files that will not be there when someone
#   configures and builds from scratch).
#
# For more information on these arguments, see `TEST_<idx>
# COPY_FILES_TO_TEST_DIR Test Blocks and Arguments
# (TRIBITS_ADD_ADVANCED_TEST())`_.
# 
# By default, each and every atomic ``TEST_<idx>`` block needs to pass (as
# defined in `Test case Pass/Fail (TRIBITS_ADD_ADVANCED_TEST())`_) in order
# for the overall test to pass.
#
# Finally, the test is only added if tests are enabled for the SE package
# (i.e. `${PACKAGE_NAME}_ENABLE_TESTS`_ ``= ON``) or the parent package (if
# this is a subpackage) (i.e. ``${PARENT_PACKAGE_NAME}_ENABLE_TESTS=ON``) or
# if other criteria are met (see some of the arguments in `Overall Arguments
# (TRIBITS_ADD_ADVANCED_TEST())`_ that can trigger a test to not be added).
# (NOTE: A more efficient way to optionally enable tests is to put them in a
# ``test/`` subdir and then include that subdir with
# `TRIBITS_ADD_TEST_DIRECTORIES()`_.)
#
# *Sections:*
#
# * `Overall Arguments (TRIBITS_ADD_ADVANCED_TEST())`_
# * `TEST_<idx> EXEC/CMND Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST())`_
# * `TEST_<idx> COPY_FILES_TO_TEST_DIR Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Test case Pass/Fail (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Overall Pass/Fail (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Argument Parsing and Ordering (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Implementation Details (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Setting Additional Test Properties (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Running multiple tests at the same time (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Disabling Tests Externally (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Debugging and Examining Test Generation (TRIBITS_ADD_ADVANCED_TEST())`_
# * `Using TRIBITS_ADD_ADVANCED_TEST() in non-TriBITS CMake projects`_
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
#   ``<testNameBase>``
#
#     The base name of the test (which will have ``${PACKAGE_NAME}_``
#     prepended to the name, see <testName> below) that will be used to name
#     the output CMake script file as well as the CTest test name passed into
#     ``ADD_TEST()``.  This must be the first argument to this function.  The
#     name is allowed to contain '/' chars but these will be replaced with
#     '__' in the overall working directory name and the ctest -P script
#     (`Debugging and Examining Test Generation
#     (TRIBITS_ADD_ADVANCED_TEST())`_).
#
#   ``OVERALL_WORKING_DIRECTORY <overallWorkingDir>``
#
#     If specified, then the working directory ``<overallWorkingDir>``
#     (relative or absolute path) will be created and all of the test commands
#     by default will be run from within this directory.  If the value
#     ``<overallWorkingDir>=TEST_NAME`` is given, then the working directory
#     will be given the name ``<testName>`` where any '/' chars are replaced
#     with '__'.  By default, if the directory ``<overallWorkingDir>`` exists
#     before the test runs, it will be deleted and created again.  If one
#     wants to preserve the contents of this directory between test runs then
#     set ``SKIP_CLEAN_OVERALL_WORKING_DIRECTORY``.  Using a separate test
#     directory is a good option to use if the commands create intermediate
#     files and one wants to make sure they get deleted before the test cases
#     are run again.  It is also important to create a separate test directory
#     if multiple tests are defined in the same ``CMakeLists.txt`` file that
#     read/write files with the same name.
#
#   ``SKIP_CLEAN_OVERALL_WORKING_DIRECTORY``
#
#     If specified, then ``<overallWorkingDir>`` will **not** be deleted if it
#     already exists.
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
#     executable command runs on.  If ``<overallNumProcs>`` is greater than
#     ``${MPI_EXEC_MAX_NUMPROCS}`` then the test will be excluded.  If not
#     specified, then the default number of processes for an MPI build will be
#     ``${MPI_EXEC_DEFAULT_NUMPROCS}``.  For serial builds, this argument is
#     ignored.  For MPI builds with all ``TEST_<idx> CMND`` blocks,
#     ``<overallNumProcs>`` is used to set the property ``PROCESSORS``. (see
#     `Running multiple tests at the same time
#     (TRIBITS_ADD_ADVANCED_TEST())`_).  **WARNING!** If just running a serial
#     script or other command, then the property ``PROCESSORS`` will still get
#     set to ``${OVERALL_NUM_MPI_PROCS}`` so in order to avoid CTest
#     unnecessarily reserving ``${OVERALL_NUM_MPI_PROCS}`` processes for a
#     serial non-MPI test, then one must leave off ``OVERALL_NUM_MPI_PROCS``
#     or explicitly pass in ``MPI_EXEC_DEFAULT_NUMPROCS 1``!
#
#   ``OVERALL_NUM_TOTAL_CORES_USED <overallNumTotalCoresUsed>``
#
#     Used for ``NUM_TOTAL_CORES_USED`` if missing in a ``TEST_<idx>`` block.
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
#   ``EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...``
#
#     If specified, gives the names of CMake variables that must evaluate to
#     true, or the test will not be added (see `TRIBITS_ADD_TEST()`_).
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
#     to run before being timed-out and killed (see `Setting timeouts for
#     tests (TRIBITS_ADD_TEST())`_).  This is for the full CTest test, not
#     individual ``TEST_<idx>`` commands!
#
#   ``ADDED_TEST_NAME_OUT <testName>``
#
#     If specified, then on output the variable ``<testName>`` will be set
#     with the name of the test passed to ``ADD_TEST()``.  Having this name
#     allows the calling ``CMakeLists.txt`` file access and set additional
#     test properties (see `Setting additional test properties
#     (TRIBITS_ADD_ADVANCED_TEST())`_).
#
# .. _TEST_<idx> EXEC/CMND Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST()):
#
# **TEST_<idx> EXEC/CMND Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST())**
#
# Each test general command block ``TEST_<idx>`` runs either a package-built
# test executable or some general command executable and is defined as either
# ``EXEC <exeRootName>`` or an arbitrary command ``CMND <cmndExec>`` with the
# arguments:
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
#     (or ``OVERALL_NUM_MPI_PROCS <overallNumProcs>`` if ``NUM_MPI_PROCS`` is
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
#     instead.  **WARNING:** If you want to run such tests using valgrind, you
#     have to use the raw executable as the ``<cmndExec>`` argument and *not*
#     the script.  For example, if you have a python script
#     ``my_python_test.py`` with ``/usr/bin/env pyhton`` at the top, you can't
#     just use::
#
#       CMND <path>/my_python_test.py ARGS <arg0> <arg1> ...
#
#     The same goes for Perl or any other scripting language.
#
#     Instead, you have to use::
#
#       CMND ${PYTHON_EXECUTABLE} ARGS <path>/my_python_test.py <arg0> <arg1> ...
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
#     If specified, then the working directory ``<workingDir>`` (relative or
#     absolute) will be created and the test will be run from within this
#     directory.  If the directory ``<workingDir>`` exists before the test
#     runs, it will be deleted and created again.  If one wants to preserve
#     the contents of this directory between test blocks, then one needs to
#     set ``SKIP_CLEAN_WORKING_DIRECTORY``.  Using a different
#     ``WORKING_DIRECTORY`` for individual test commands allows creating
#     independent working directories for each test case.  This would be
#     useful if a single ``OVERALL_WORKING_DIRECTORY`` was not sufficient for
#     some reason.
#
#   ``SKIP_CLEAN_WORKING_DIRECTORY``
#
#     If specified, then ``<workingDir>`` will **not** be deleted if it
#     already exists.
#
#   ``NUM_MPI_PROCS <numProcs>``
#
#     If specified, then ``<numProcs>`` is the number of processors used for
#     MPI executables.  If not specified, this will default to
#     ``<overallNumProcs>`` from ``OVERALL_NUM_MPI_PROCS <overallNumProcs>``.
#
#   ``NUM_TOTAL_CORES_USED <numTotalCoresUsed>``
#
#     If specified, gives the total number of processes used by this
#     command/executable.  If this is missing, but ``NUM_MPI_PROCS
#     <numProcs>`` is specified, then ``<numProcs>`` is used instead.  If
#     ``NUM_TOTAL_CORES_USED`` is missing BUT ``OVERALL_NUM_TOTAL_CORES_USED
#     <overallNumTotalCoresUsed>`` is, then ``<overallNumTotalCoresUsed>`` is
#     used for ``<numTotalCoresUsed>``.  This argument is used for test
#     scripts/executables that use more cores than MPI processes
#     (i.e. ``<numProcs>``) and its only purpose is to inform CTest and
#     TriBITS of the maximum number of cores that are used by the underlying
#     test executable/script.  When ``<numTotalCoresUsed>`` is greater than
#     ``${MPI_EXEC_MAX_NUMPROCS}``, then the test will not be added.
#     Otherwise, the CTest property ``PROCESSORS`` is set to the max over all
#     ``<numTotalCoresUsed>`` so that CTest knows how to best schedule the
#     test w.r.t. other tests on a given number of available processes.
#
#   ``OUTPUT_FILE <outputFile>``
#
#     If specified, then stdout and stderr for the test case will be sent to
#     ``<outputFile>``.  By default, the contents of this file will **also**
#     be printed to STDOUT unless ``NO_ECHO_OUT`` is passed as well.
#
#     NOTE: Contrary to CMake documentation for EXECUTE_PROCESS(), STDOUT and
#     STDERR may not get output in the correct order interleaved correctly,
#     even in serial without MPI.  Therefore, you can't write any tests that
#     depend on the order of STDOUT and STDERR output in relation to each
#     other.  Also note that all of STDOUT and STDERR will be first read into
#     the CTest executable process main memory before the file
#     ``<outputFile>`` is written.  Therefore, don't run executables or
#     commands that generate massive amounts of console output or it may
#     exhaust main memory.  Instead, have the command or executable write
#     directly to a file instead of going through STDOUT.
#
#   ``NO_ECHO_OUTPUT``
#
#     If specified, then the output for the test command will not be echoed to
#     the output for the entire test command.
#
# By default, an individual test case ``TEST_<IDX>`` is assumed to pass if the
# executable or commands returns a non-zero value to the shell.  However, a
# test case can also be defined to pass or fail based on the arguments/options
# (see `Test case Pass/Fail (TRIBITS_ADD_ADVANCED_TEST())`_):
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
#     given regular expression.  Otherwise, it is assumed to fail.  TIPS:
#     Replace ';' with '[;]' or CMake will interpret this as an array elemnt
#     boundary.  To match '.', use '[.]'.
#
#   ``PASS_REGULAR_EXPRESSION_ALL "<regex1>" "<regex2>" ... "<regexn>"``
#
#     If specified, the test command will be assumed to pass if the output
#     matches all of the provided regular expressions.  Note that this is not
#     a capability of raw ctest and represents an extension provided by
#     TriBITS.  NOTE: It is critical that you replace ';' with '[;]' or CMake
#     will interpret this as an array elemnt boundary.
#
#   ``STANDARD_PASS_OUTPUT``
#
#     If specified, the test command will be assumed to pass if the string
#     expression "Final Result: PASSED" is found in the output for the test.
#     This as the result of directly passing in ``PASS_REGULAR_EXPRESSION
#     "Final Result: PASSED"``.
#
#   ``FAIL_REGULAR_EXPRESSION "<regex>"``
#
#     If specified, the test command will be assumed to fail if it matches the
#     given regular expression.  Otherwise, it is assumed to pass.  This will
#     be applied and take precedence over other above pass criteria.  For
#     example, if even if ``PASS_REGULAR_EXPRESSION`` or
#     ``PASS_REGULAR_EXPRESSION_ALL`` match, then the test will be marked as
#     failed if this fail regex matches the output.
#
#   ``ALWAYS_FAIL_ON_NONZERO_RETURN``
#
#     If specified, then the test case will be marked as failed if the test
#     command returns nonzero, independent of the other pass/fail criteria.
#     This option is used in cases where one wants to grep for strings in the
#     output but still wants to require a zero return code.  This make for a
#     stronger test by requiring that both the strings are found and that the
#     command returns 0.
#
#   ``ALWAYS_FAIL_ON_ZERO_RETURN``
#
#     If specified, then the test case will be marked as failed if the test
#     command returns zero '0', independent of the other pass/fail criteria.
#     This option is used in cases where one wants to grep for strings in the
#     output but still wants to require a nonzero return code.  This make for
#     a stronger test by requiring that both the strings are found and that
#     the command returns != 0.
#
#   ``WILL_FAIL``
#
#      If specified, invert the result from the other pass/fail criteria.  For
#      example, if the regexes in ``PASS_REGULAR_EXPRESSION`` or
#      ``PASS_REGULAR_EXPRESSION_ALL`` indicate that a test should pass, then
#      setting ``WILL_FAIL`` will invert that and report the test as failing.
#      But typically this is used to report a test that returns a nonzero code
#      as passing.
# 
# All of the arguments for a test block ``TEST_<idx>`` must appear directly
# below their ``TEST_<idx>`` argument and before the next test block (see
# `Argument Parsing and Ordering (TRIBITS_ADD_ADVANCED_TEST())`_).
#
# **WARNING:** The current implementation limits the number of ``TEST_<idx>``
# cases to just 20 (i.e. ``<idx>=0...19``).  And if more test cases are added
# (e.g. ``TEST_20``), the current implementation can't detect that case and
# the resulting behavior is undefined.  This restriction will be removed in a
# future version of TriBITS.
#
# .. _TEST_<idx> COPY_FILES_TO_TEST_DIR Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST()):
#
# **TEST_<idx> COPY_FILES_TO_TEST_DIR Test Blocks and Arguments (TRIBITS_ADD_ADVANCED_TEST())**
#
# The arguments for the ``TEST_<idx>`` ``COPY_FILES_TO_TEST_DIR`` block are:
#
#   ``COPY_FILES_TO_TEST_DIR <file0> <file1> ... <filen>``
#
#     Required list of 1 or more file names for files that will be copied from
#     ``<srcDir>/`` to ``<destDir>/``.
#
#   ``SOURCE_DIR <srcDir>``
#
#     Optional source directory where the files will be copied from.  If
#     ``<srcDir>`` is not given, then it is assumed to be
#     ``${CMAKE_CURRENT_SOURCE_DIR}``.  If ``<srcDir>`` is given but is a
#     relative path, then it is interpreted relative to
#     ``${CMAKE_CURRENT_SOURCE_DIR}``.  If ``<srcDir>`` is an absolute path,
#     then that path is used without modification.
#
#   ``DEST_DIR <destDir>``
#
#     Optional destination directory where the files will be copied to.  If
#     ``<destDir>`` is not given, then it is assumed to be the working
#     directory where the test is running (typically a new directory created
#     under ``${CMAKE_CURRENT_BINARY_DIR}`` when ``OVERALL_WORKING_DIR
#     TEST_NAME`` is given).  If ``<destDir>`` is given but is a relative
#     path, then it is interpreted relative to the current test working
#     directory.  If ``<destDir>`` is an absolute path, then that path is used
#     without modification.  If ``<destDir>`` does not exist, then it will be
#     created (including several directory levels deep if needed).
#
# .. _Test case Pass/Fail (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Test case Pass/Fail (TRIBITS_ADD_ADVANCED_TEST())**
#
# The logic for how pass/fail for a ``TEST_<IDX>`` ``EXEC`` or ``CMND`` case
# is applied is given by::
#
#   # A) Apply first set of pass/fail logic
#   TEST_CASE_PASSED = FALSE
#   If PASS_ANY specified:
#     TEST_CASE_PASSED = TRUE
#   Else If PASS_REGULAR_EXPRESSION specified and "<regex>" matches:
#     TEST_CASE_PASSED = TRUE
#   Else if PASS_REGULAR_EXPRESSION_ALL specified:
#     TEST_CASE_PASSED = TRUE
#     For each "<regexi>":
#       If "<regixi>" does not match:
#         TEST_CASE_PASSED = FALSE
#   Else
#     If command return code == 0:
#       TEST_CASE_PASSED = TRUE
#     Endif
#   Endif
#
#   # B) Check for failing regex matching?
#   If FAIL_REGULAR_EXPRESSION specified and "<regex>" matches:
#     TEST_CASE_PASSED = FALSE
#   Endif
#
#   # C) Check for return code always 0 or !=0?
#   If ALWAYS_FAIL_ON_NONZERO_RETURN specified and return code != 0:
#     TEST_CASE_PASSED = FALSE
#   ElseIf ALWAYS_FAIL_ON_ZERO_RETURN specified and return code == 0:
#     TEST_CASE_PASSED = FALSE
#   Endif
#
#   # D) Invert pass/fail result?
#   If WILL_FAIL specified:
#     If TEST_CASE_PASSED:
#       TEST_CASE_PASSED = FALSE
#     Else
#       TEST_CASE_PASSED = TRUE
#     Endif
#   Endif
#
# .. _Overall Pass/Fail (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Overall Pass/Fail (TRIBITS_ADD_ADVANCED_TEST())**
#
# By default, the overall test will be assumed to pass if it prints::
#
#   "OVERALL FINAL RESULT: TEST PASSED (<testName>)"
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
# ``CMAKE_PARSE_ARGUMENTS()`` which has a certain set of behaviors.  The
# parsing using ``CMAKE_PARSE_ARGUMENTS()`` is actually done in two phases.
# There is a top-level parsing of the "overall" arguments listed in `Overall
# Arguments (TRIBITS_ADD_ADVANCED_TEST())`_ that also pulls out the test
# blocks.  Then there is a second level of parsing using
# ``CMAKE_PARSE_ARGUMENTS()`` for each of the ``TEST_<idx>`` blocks.  Because
# of this usage, there are a few restrictions that one needs to be aware of
# when using ``TRIBITS_ADD_ADVANCED_TEST()``.  This short sections tries to
# explain the behaviors and what is allowed and what is not allowed.
#
# For the most part, the "overall" arguments and the arguments inside of any
# individual ``TEST_<idx>`` blocks can be listed in any order but there are
# restrictions related to the grouping of overall arguments and ``TEST_<idx>``
# blocks which are as follows:
#
# * The ``<testNameBase>`` argument must be the first listed (it is the only
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
# name ``<testName>.cmake`` (with any '/' replaced with '__') gets created in
# the current binary directory that then gets added to CTest using::
#
#   ADD_TEST(<testName> cmake [other options] -P <testName>.cmake)
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
# ``SET_TESTS_PROPERTIES(<testName> ...)``, where ``<testName>`` is returned
# using the ``ADDED_TEST_NAME_OUT <testName>`` argument.  Therefore, any tests
# properties that are not directly supported by this function and passed
# through the argument list to this wrapper function can be set in the outer
# ``CMakeLists.txt`` file after the call to ``TRIBITS_ADD_ADVANCED_TEST()``.
# For example::
#
#   TRIBITS_ADD_ADVANCED_TEST_TEST( someTest ...
#     ADDED_TEST_NAME_OUT  someTest_TEST_NAME )
#
#   IF (someTest_TEST_NAME)
#     SET_TESTS_PROPERTIES( ${someTest_TEST_NAME}
#       PROPERTIES ATTACHED_FILES someTest.log )
#   ENDIF()
#
# where the test writes a log file ``someTest.log`` that we want to submit to
# CDash also.
#
# This approach will work no matter what TriBITS names the individual test(s)
# or whether the test(s) are added or not (depending on other arguments like
# ``COMM``, ``XHOST``, etc.).
#
# The following built-in CTest test properties are set through `Overall
# Arguments (TRIBITS_ADD_ADVANCED_TEST())`_ or are otherwise automatically set
# by this function and should **NOT** be overridden by direct calls to
# ``SET_TESTS_PROPERTIES()``: ``ENVIRONMENT``, ``FAIL_REGULAR_EXPRESSION``,
# ``LABELS``, ``PASS_REGULAR_EXPRESSION``, ``RUN_SERIAL``, ``TIMEOUT``,
# ``WILL_FAIL``, and ``WORKING_DIRECTORY``.
#
# However, generally, other built-in CTest test properties can be set after
# the test is added like show above.  Examples of test properties that can be
# set using direct calls to ``SET_TESTS_PROPERTIES()`` include
# ``ATTACHED_FILES``, ``ATTACHED_FILES_ON_FAIL``, ``COST``, ``DEPENDS``,
# ``MEASUREMENT``, and ``RESOURCE_LOCK``.
#
# For example, one can set a dependency between two tests using::
#
#   TRIBITS_ADD_ADVANCED_TEST_TEST( test_a [...]
#      ADDED_TEST_NAME_OUT  test_a_TEST_NAME )
#   
#   TRIBITS_ADD_ADVANCED_TEST_TEST( test_b [...]
#      ADDED_TEST_NAME_OUT  test_z_TEST_NAME )
#   
#   IF (test_a_TEST_NAME AND test_b_TEST_NAME)
#     SET_TESTS_PROPERTIES(${test_b_TEST_NAME}
#       PROPERTIES DEPENDS ${test_a_TEST_NAME})
#   ENDIF()
#
# This ensures that test ``test_b`` will always be run after ``test_a`` if
# both tests are run by CTest.
#
# .. _Running multiple tests at the same time (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Running multiple tests at the same time (TRIBITS_ADD_ADVANCED_TEST())**
#
# Just as with `TRIBITS_ADD_TEST()`_, setting ``NUM_MPI_PROCS <numProcs>`` or
# ``OVERALL_NUM_MPI_PROCS <numOverallProcs>`` or ``NUM_TOTAL_CORES_USED
# <numTotalCoresUsed>`` or ``OVERALL_NUM_TOTAL_CORES_USED
# <overallNumTotalCoresUsed>`` will set the ``PROCESSORS`` CTest property to
# allow CTest to schedule and run multiple tests at the same time when ``'ctest
# -j<N>'`` is used (see `Running multiple tests at the same time
# (TRIBITS_ADD_TEST())`_).
#
# .. _Disabling Tests Externally (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Disabling Tests Externally (TRIBITS_ADD_ADVANCED_TEST())**
#
# The test can be disabled externally by setting the CMake cache variable
# ``<testName>_DISABLE=TRUE``.  This allows tests to be disabled on a
# case-by-case basis.  The name ``<testName>`` must be the *exact* name
# that shows up in ``ctest -N`` when running the test.
#
# .. _Debugging and Examining Test Generation (TRIBITS_ADD_ADVANCED_TEST()):
#
# **Debugging and Examining Test Generation (TRIBITS_ADD_ADVANCED_TEST())**
#
# In order to see what tests get added and if not then why, configure with
# ``${PROJECT_NAME}_TRACE_ADD_TEST=ON``.  That will print one line per test
# that shows that the test got added or not and if not then why the test was
# not added (i.e. due to ``COMM``, ``OVERALL_NUM_MPI_PROCS``,
# ``NUM_MPI_PROCS``, ``CATEGORIES``, ``HOST``, ``XHOST``, ``HOSTTYPE``, or
# ``XHOSTTYPE``).
#
# Likely the best way to debug test generation using this function is to
# examine the generated file ``<testName>.cmake`` in the current binary
# directory (see `Implementation Details (TRIBITS_ADD_ADVANCED_TEST())`_) and
# the generated ``CTestTestfile.cmake`` file that should list this test case.
#
# .. _Using TRIBITS_ADD_ADVANCED_TEST() in non-TriBITS CMake projects:
#
# **Using TRIBITS_ADD_ADVANCED_TEST() in non-TriBITS CMake projects**
#
# The function ``TRIBITS_ADD_ADVANCED_TEST()`` can be used to add tests in
# non-TriBITS projects.  To do so, one just needs to set the variables
# ``PROJECT_NAME``, ``PACKAGE_NAME`` (which could be the same as
# ``PROJECT_NAME``), ``${PACKAGE_NAME}_ENABLE_TESTS=TRUE``, and
# ``${PROJECT_NAME}_TRIBITS_DIR`` (pointing to the TriBITS location).  For example,
# a valid project can be a simple as::
#
#   CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
#   SET(PROJECT_NAME TAATDriver)
#   PROJECT(${PROJECT_NAME} NONE)
#   SET(${PROJECT_NAME}_TRACE_ADD_TEST TRUE)
#   SET(${PROJECT_NAME}_TRIBITS_DIR ""  CACHE FILEPATH
#     "Location of TriBITS to use." ) 
#   SET(PACKAGE_NAME ${PROJECT_NAME})
#   SET(${PACKAGE_NAME}_ENABLE_TESTS TRUE)
#   SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
#     ${TRIBITS_DIR}/core/utils
#     ${TRIBITS_DIR}/core/package_arch )
#   INCLUDE(TribitsAddAdvancedTest)
#   INCLUDE(CTest)
#   ENABLE_TESTING()
#   
#   TRIBITS_ADD_ADVANCED_TEST(
#     TAAT_COPY_FILES_TO_TEST_DIR_bad_file_name
#     OVERALL_WORKING_DIRECTORY TEST_NAME
#     TEST_0 CMND echo ARGS "Hello World!"
#       PASS_REGULAR_EXPRESIOIN "Hello World"
#     )
#
FUNCTION(TRIBITS_ADD_ADVANCED_TEST TEST_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_ADD_ADVANCED_TEST: ${TEST_NAME_IN}\n")
  ENDIF()

  GLOBAL_SET(TRIBITS_SET_TEST_PROPERTIES_INPUT)
  GLOBAL_SET(MESSAGE_WRAPPER_INPUT)

  # Set the full TEST_NAME
  IF (PACKAGE_NAME)
    SET(TEST_NAME ${PACKAGE_NAME}_${TEST_NAME_IN})
  ELSE()
    SET(TEST_NAME ${TEST_NAME_IN})
  ENDIF()

  IF (NOT CMAKE_VERSION VERSION_LESS "3.1.0.0")
    # Avoid quoted strings lookup variables 
    CMAKE_POLICY(SET CMP0054 NEW)
    # NOTE: For some reason, setting this policy at the top level with TriBITS
    # in TribitsCMakePolices.cmake does not affect this function.  Therefore,
    # I have to set it again here.
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

  CMAKE_PARSE_ARGUMENTS(
     #prefix
     PARSE
     #options
     "FAIL_FAST;RUN_SERIAL;SKIP_CLEAN_OVERALL_WORKING_DIRECTORY"
     # one_value_keywords
     ""
     # multi_value_keywords
     "${TEST_IDX_LIST};OVERALL_WORKING_DIRECTORY;KEYWORDS;COMM;OVERALL_NUM_MPI_PROCS;OVERALL_NUM_TOTAL_CORES_USED;FINAL_PASS_REGULAR_EXPRESSION;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;EXCLUDE_IF_NOT_TRUE;FINAL_FAIL_REGULAR_EXPRESSION;TIMEOUT;ENVIRONMENT;ADDED_TEST_NAME_OUT"
     ${ARGN}
     )

  IF(PARSE_ADDED_TEST_NAME_OUT)
    SET(${PARSE_ADDED_TEST_NAME_OUT} "" PARENT_SCOPE )
  ENDIF()

  # Set the name of the cmake -P script.
  STRING(REPLACE "/" "__" NORMALIZED_TEST_NAME "${TEST_NAME}")
  SET(TEST_SCRIPT_FILE_NAME "${NORMALIZED_TEST_NAME}.cmake")

  # Set the relative overall working directory and abs working directory
  IF (PARSE_OVERALL_WORKING_DIRECTORY)
    IF ("${PARSE_OVERALL_WORKING_DIRECTORY}" STREQUAL "TEST_NAME")
      SET(PARSE_OVERALL_WORKING_DIRECTORY ${NORMALIZED_TEST_NAME})
    ENDIF()
   # Test will run in created working subdir
   SET(ABS_OVERALL_WORKING_DIRECTORY
     ${CMAKE_CURRENT_BINARY_DIR}/${PARSE_OVERALL_WORKING_DIRECTORY})
  ELSE()
    # Test runs in current binary directory (not a good idea!) 
    SET(ABS_OVERALL_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()

  #
  # B) Add or don't add tests based on a number of criteria
  #

  SET(ADD_THE_TEST FALSE)
  TRIBITS_ADD_TEST_PROCESS_ENABLE_TESTS(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

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

  IF (PARSE_OVERALL_NUM_MPI_PROCS  AND  PARSE_OVERALL_NUM_TOTAL_CORES_USED)
    IF (PARSE_OVERALL_NUM_MPI_PROCS  GREATER  PARSE_OVERALL_NUM_TOTAL_CORES_USED)
      MESSAGE_WRAPPER(FATAL_ERROR
        "ERROR: ${TEST_NAME}: OVERALL_NUM_MPI_PROCS='${PARSE_OVERALL_NUM_MPI_PROCS}' > OVERALL_NUM_TOTAL_CORES_USED='${PARSE_OVERALL_NUM_TOTAL_CORES_USED}' not allowed!")
      RETURN()
    ENDIF()
  ENDIF()

  # ToDo: Assert that OVERALL_NUM_TOTAL_CORES_USED >= OVERALL_NUM_MPI_PROCS

  IF (PARSE_OVERALL_NUM_MPI_PROCS)
    SET(MAX_NUM_MPI_PROCS_USED ${PARSE_OVERALL_NUM_MPI_PROCS})
    SET(MAX_NUM_PROCESSORS_USED ${PARSE_OVERALL_NUM_MPI_PROCS})
  ELSE()
    SET(MAX_NUM_MPI_PROCS_USED ${PARSE_OVERALL_NUM_MPI_PROCS})
    SET(MAX_NUM_PROCESSORS_USED 1)
  ENDIF()

  SET(HAS_AT_LEAST_ONE_EXEC FALSE)

  FOREACH( TEST_CMND_IDX RANGE ${MAX_NUM_TEST_CMND_IDX} )

    IF (NOT PARSE_TEST_${TEST_CMND_IDX} )
      BREAK()
    ENDIF()

    IF (NOT  ADD_THE_TEST)
      BREAK()
    ENDIF()

    MATH( EXPR NUM_CMNDS ${NUM_CMNDS}+1 )

    # Parse the test command case

    #PRINT_VAR(PARSE_TEST_${TEST_CMND_IDX})

    # Search to see if we are copying files or not for this TEST_<IDX> block ...

    SET(PARSE_COPY_FILES_TO_TEST_DIR)
    SET(COPY_FILES_TO_TEST_DIR_IDX FALSE)
    FOREACH(PARSE_TEST_IDX_ARGS ${PARSE_TEST_${TEST_CMND_IDX}})
      IF (PARSE_TEST_IDX_ARGS STREQUAL "COPY_FILES_TO_TEST_DIR")
        SET(COPY_FILES_TO_TEST_DIR_IDX TRUE)
      ENDIF()
    ENDFOREACH()

    IF (COPY_FILES_TO_TEST_DIR_IDX)

      # Do a special parse just for TEST_<IDX> blocks of type
      # COPY_FILES_TO_TEST_DIR

      CMAKE_PARSE_ARGUMENTS(
         #prefix
         PARSE
         #options
         ""
         # one_value_keywords
         ""
         # multi_value_keywords
         "COPY_FILES_TO_TEST_DIR;SOURCE_DIR;DEST_DIR" 
	 # Arguments to parse
         ${PARSE_TEST_${TEST_CMND_IDX}}
         )
      TRIBITS_CHECK_FOR_UNPARSED_ARGUMENTS()
      TRIBITS_ASSERT_PARSE_ARG_ONE_OR_MORE_VALUES(PARSE COPY_FILES_TO_TEST_DIR)
      TRIBITS_ASSERT_PARSE_ARG_ZERO_OR_ONE_VALUE(PARSE SOURCE_DIR)
      TRIBITS_ASSERT_PARSE_ARG_ZERO_OR_ONE_VALUE(PARSE DEST_DIR)

      SET(PARSE_EXEC)
      SET(PARSE_CMND)

    ELSE()

      # Parse TEST_<IDX> block args for types EXEC and CMND

      CMAKE_PARSE_ARGUMENTS(
         #prefix
         PARSE
         #options
          "NOEXEPREFIX;NOEXESUFFIX;NO_ECHO_OUTPUT;PASS_ANY;STANDARD_PASS_OUTPUT;ALWAYS_FAIL_ON_NONZERO_RETURN;ALWAYS_FAIL_ON_ZERO_RETURN;WILL_FAIL;ADD_DIR_TO_NAME;SKIP_CLEAN_WORKING_DIRECTORY"
         # one_value_keywords
         ""
         # multi_value_keywords
         "EXEC;CMND;ARGS;DIRECTORY;MESSAGE;WORKING_DIRECTORY;OUTPUT_FILE;NUM_MPI_PROCS;NUM_TOTAL_CORES_USED;PASS_REGULAR_EXPRESSION_ALL;FAIL_REGULAR_EXPRESSION;PASS_REGULAR_EXPRESSION"
         ${PARSE_TEST_${TEST_CMND_IDX}}
         )
  
      TRIBITS_CHECK_FOR_UNPARSED_ARGUMENTS()

    ENDIF()

    #
    # Set up the command that will be written into the cmake -P *.cmake file
    #

    SET(ARGS_STR ${PARSE_ARGS})
    #PRINT_VAR(ARGS_STR)
    #IF (PARSE_ARGS)
    #  TRIBITS_JOIN_EXEC_PROCESS_SET_ARGS( ARGS_STR ${PARSE_ARGS} )
    #ENDIF()

    IF (PARSE_EXEC)

      #
      # This is an EXEC test block
      #

      SET(HAS_AT_LEAST_ONE_EXEC TRUE)

      LIST( LENGTH PARSE_EXEC PARSE_EXEC_LEN )
      IF (NOT PARSE_EXEC_LEN EQUAL 1)
        MESSAGE(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} EXEC = '${PARSE_EXEC}'"
          " must be a single name.  To add arguments use ARGS <arg1> <arg2> ...." )
      ENDIF()

      TRIBITS_ADD_TEST_GET_EXE_BINARY_NAME( "${PARSE_EXEC}"
        ${PARSE_NOEXEPREFIX} ${PARSE_NOEXESUFFIX}
        ${PARSE_ADD_DIR_TO_NAME} EXE_BINARY_NAME )

      TRIBITS_ADD_TEST_ADJUST_DIRECTORY( ${EXE_BINARY_NAME} "${PARSE_DIRECTORY}"
        EXECUTABLE_PATH)

      IF (PARSE_NUM_MPI_PROCS)
        SET(NUM_MPI_PROC_VAR_NAME "NUM_MPI_PROCS")
      ELSE()
        SET(PARSE_NUM_MPI_PROCS ${PARSE_OVERALL_NUM_MPI_PROCS})
        SET(NUM_MPI_PROC_VAR_NAME "OVERALL_NUM_MPI_PROCS")
      ENDIF()

      TRIBITS_ADD_TEST_GET_NUM_PROCS_USED("${PARSE_NUM_MPI_PROCS}"
        "${NUM_MPI_PROC_VAR_NAME}"  NUM_PROCS_USED  NUM_PROCS_USED_NAME)

      IF (NUM_PROCS_USED LESS 0)
        SET(ADD_THE_TEST FALSE)
      ELSEIF (NUM_PROCS_USED  GREATER  MAX_NUM_MPI_PROCS_USED)
        SET(MAX_NUM_MPI_PROCS_USED  ${NUM_PROCS_USED})
      ENDIF()

      IF (PARSE_NUM_TOTAL_CORES_USED)
        SET(NUM_TOTAL_CORES_USED  ${PARSE_NUM_TOTAL_CORES_USED})
        SET(NUM_TOTAL_CORES_USED_NAME  "NUM_TOTAL_CORES_USED")
      ELSE()
        SET(NUM_TOTAL_CORES_USED  ${PARSE_OVERALL_NUM_TOTAL_CORES_USED})
        SET(NUM_TOTAL_CORES_USED_NAME  "OVERALL_NUM_TOTAL_CORES_USED")
      ENDIF()

      TRIBITS_ADD_TEST_GET_NUM_TOTAL_CORES_USED("${TEST_NAME}${MPI_NAME_POSTFIX}"
        "${NUM_TOTAL_CORES_USED}"  "${NUM_TOTAL_CORES_USED_NAME}"
        "${NUM_PROCS_USED}"  "${NUM_PROCS_USED_NAME}"
        NUM_TOTAL_CORES_USED  SKIP_TEST)
      IF (SKIP_TEST)
        SET(ADD_THE_TEST FALSE)
      ENDIF()

      IF (NUM_TOTAL_CORES_USED  GREATER  MAX_NUM_PROCESSORS_USED)
        SET(MAX_NUM_PROCESSORS_USED  ${NUM_TOTAL_CORES_USED})
      ENDIF()

      IF(ADD_THE_TEST)
        LIST(APPEND TEST_EXE_LIST ${EXECUTABLE_PATH})
      ENDIF()

      TRIBITS_ADD_TEST_GET_TEST_CMND_ARRAY( TEST_CMND_ARRAY
        "${EXECUTABLE_PATH}" "${NUM_PROCS_USED}" ${ARGS_STR} )
      #PRINT_VAR(TEST_CMND_ARRAY)

    ELSEIF (PARSE_CMND)

      #
      # This is a COMMAND test block
      #

      LIST( LENGTH PARSE_CMND PARSE_CMND_LEN )
      IF (NOT PARSE_CMND_LEN EQUAL 1)
        MESSAGE(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} CMND = '${PARSE_CMND}'"
          " must be a single command.  To add arguments use ARGS <arg1> <arg2> ...." )
      ENDIF()

      IF (PARSE_NUM_TOTAL_CORES_USED)
        SET(NUM_TOTAL_CORES_USED  ${PARSE_NUM_TOTAL_CORES_USED})
        SET(NUM_TOTAL_CORES_USED_NAME  "NUM_TOTAL_CORES_USED")
      ELSE()
        SET(NUM_TOTAL_CORES_USED  ${PARSE_OVERALL_NUM_TOTAL_CORES_USED})
        SET(NUM_TOTAL_CORES_USED_NAME  "OVERALL_NUM_TOTAL_CORES_USED")
      ENDIF()

      TRIBITS_ADD_TEST_GET_NUM_TOTAL_CORES_USED("${TEST_NAME}${MPI_NAME_POSTFIX}"
        "${NUM_TOTAL_CORES_USED}"  "${NUM_TOTAL_CORES_USED_NAME}"
        "1"  "DUMMY_NUM_MPI_PROCS"  # Never be printed
        NUM_TOTAL_CORES_USED  SKIP_TEST)
      IF (SKIP_TEST)
        SET(ADD_THE_TEST FALSE)
      ENDIF()

      IF (NUM_TOTAL_CORES_USED  GREATER  MAX_NUM_PROCESSORS_USED)
        SET(MAX_NUM_PROCESSORS_USED  ${NUM_TOTAL_CORES_USED})
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

    ELSEIF (PARSE_COPY_FILES_TO_TEST_DIR)

      #
      # This is a COPY_FLES_TO_TEST_DIR block
      #

      # FILES_TO_COPY_COMMA_SEP
      SET(FILES_TO_COPY_COMMA_SEP "${PARSE_COPY_FILES_TO_TEST_DIR}")
      string(REPLACE ";" "," FILES_TO_COPY_COMMA_SEP
        "${FILES_TO_COPY_COMMA_SEP}" )
      # NOTE: Above, we have to replace ';' with ',' or the lower commands
      # APPEND_STRING_VAR() will replace ';' with ''.  This is *not* what we
      # want.  In DriveAdvancedTest.cmake, we will replace the ',' with ';'
      # again :-)  

      # SOURCE_DIR
      IF (PARSE_SOURCE_DIR)
        IF (IS_ABSOLUTE "${PARSE_SOURCE_DIR}")
          SET(COPY_FILES_TO_TEST_DIR_SOURCE_DIR
            "${PARSE_SOURCE_DIR}")
        ELSE()
          SET(COPY_FILES_TO_TEST_DIR_SOURCE_DIR
            "${CMAKE_CURRENT_SOURCE_DIR}/${PARSE_SOURCE_DIR}")
        ENDIF()
      ELSE()
        SET(COPY_FILES_TO_TEST_DIR_SOURCE_DIR
          "${CMAKE_CURRENT_SOURCE_DIR}")
      ENDIF()

      # DEST_DIR
      IF (PARSE_DEST_DIR)
        IF (IS_ABSOLUTE "${PARSE_DEST_DIR}")
          SET(COPY_FILES_TO_TEST_DIR_DEST_DIR
            "${PARSE_DEST_DIR}")
        ELSE()
          SET(COPY_FILES_TO_TEST_DIR_DEST_DIR
            "${ABS_OVERALL_WORKING_DIRECTORY}/${PARSE_DEST_DIR}")
        ENDIF()
      ELSE()
        SET(COPY_FILES_TO_TEST_DIR_DEST_DIR
          "${ABS_OVERALL_WORKING_DIRECTORY}")
      ENDIF()

    ELSE()

      MESSAGE( FATAL_ERROR
        "Must have EXEC, CMND, or COPY_FILES_TO_TEST_DIR for TEST_${TEST_CMND_IDX}" )

    ENDIF()

    #
    # Write parts for this TEST_<IDX> block to TEST_SCRIPT_STR
    #

    IF (PARSE_COPY_FILES_TO_TEST_DIR)

      # Write the vars for COPY_FILES_TO_TEST_DIR 
  
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_COPY_FILES_TO_TEST_DIR"
        " \"${FILES_TO_COPY_COMMA_SEP}\")\n"
        )
      IF (TRIBITS_ADD_ADVANCED_TEST_UNITTEST)
        GLOBAL_SET(TRIBITS_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX}
          "${TEST_CMND_STR}" )
      ENDIF()
  
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_SOURCE_DIR"
        " \"${COPY_FILES_TO_TEST_DIR_SOURCE_DIR}\")\n"
        )
  
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_DEST_DIR"
        " \"${COPY_FILES_TO_TEST_DIR_DEST_DIR}\")\n"
        )

    ELSE()

      # Write the command to be run for EXEC and CMND blocks ...

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
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_SKIP_CLEAN_WORKING_DIRECTORY ${PARSE_SKIP_CLEAN_WORKING_DIRECTORY} )\n"
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

    IF (PARSE_FAIL_REGULAR_EXPRESSION)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_FAIL_REGULAR_EXPRESSION \"${PARSE_FAIL_REGULAR_EXPRESSION}\" )\n"
        )
    ENDIF()

    IF (PARSE_ALWAYS_FAIL_ON_NONZERO_RETURN)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_ALWAYS_FAIL_ON_NONZERO_RETURN TRUE )\n"
        )
    ENDIF()

    IF (PARSE_ALWAYS_FAIL_ON_ZERO_RETURN)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_ALWAYS_FAIL_ON_ZERO_RETURN TRUE )\n"
        )
    ENDIF()

    IF (PARSE_WILL_FAIL)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_${TEST_CMND_IDX}_WILL_FAIL TRUE )\n"
        )
    ENDIF()

  ENDFOREACH()

  # ToDo: Verify that TEST_${MAX_NUM_TEST_CMND_IDX}+1 does *not* exist!

  #
  # F) Set the CTest test to run the new script
  #

  IF (ADD_THE_TEST)

    #
    # F.1) Call ADD_TEST() and set the test properties
    #
  
    SET(TEST_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_SCRIPT_FILE_NAME}")

    IF(NOT TRIBITS_ADD_TEST_ADD_TEST_UNITTEST)
      # Tell CTest to run our script for this test.  Pass the test-type
      # configuration name to the script in the TEST_CONFIG variable.
      ADD_TEST( ${TEST_NAME}
        ${CMAKE_COMMAND} "-DTEST_CONFIG=\${CTEST_CONFIGURATION_TYPE}"
        -P "${TEST_SCRIPT_FILE}")
    ENDIF()

    IF(PARSE_ADDED_TEST_NAME_OUT)
      SET(${PARSE_ADDED_TEST_NAME_OUT} ${TEST_NAME} PARENT_SCOPE )
    ENDIF()

    LIST(REMOVE_DUPLICATES TEST_EXE_LIST)
    TRIBITS_SET_TEST_PROPERTY(${TEST_NAME} PROPERTY REQUIRED_FILES ${TEST_EXE_LIST})

    IF(PARSE_RUN_SERIAL)
      TRIBITS_SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES RUN_SERIAL ON)
    ENDIF()

    TRIBITS_PRIVATE_ADD_TEST_ADD_LABEL_AND_KEYWORDS(${TEST_NAME})

    #This if clause will set the number of PROCESSORS to reserve during testing
    #to the number requested for the test.
    IF(MAX_NUM_PROCESSORS_USED)
      TRIBITS_SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES
        PROCESSORS "${MAX_NUM_PROCESSORS_USED}")
    ENDIF()

    IF (PARSE_FINAL_PASS_REGULAR_EXPRESSION)
      TRIBITS_SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        PASS_REGULAR_EXPRESSION "${PARSE_FINAL_PASS_REGULAR_EXPRESSION}" )
    ELSEIF (PARSE_FINAL_FAIL_REGULAR_EXPRESSION)
      TRIBITS_SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        FAIL_REGULAR_EXPRESSION "${PARSE_FINAL_FAIL_REGULAR_EXPRESSION}" )
    ELSE()
      TRIBITS_SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        PASS_REGULAR_EXPRESSION
        "OVERALL FINAL RESULT: TEST PASSED .${TEST_NAME}." )
    ENDIF()

    TRIBITS_PRIVATE_ADD_TEST_SET_TIMEOUT(${TEST_NAME}  TIMEOUT_USED)

    TRIBITS_PRIVATE_ADD_TEST_SET_ENVIRONMENT(${TEST_NAME})

    IF (TPL_ENABLE_MPI AND HAS_AT_LEAST_ONE_EXEC)
      SET(MAX_NUM_MPI_PROCS_USED_TO_PRINT  ${MAX_NUM_MPI_PROCS_USED})
    ELSE()
      SET(MAX_NUM_MPI_PROCS_USED_TO_PRINT "")
    ENDIF()

    TRIBITS_PRIVATE_ADD_TEST_PRINT_ADDED(${TEST_NAME}
      "${PARSE_CATEGORIES}"  "${MAX_NUM_MPI_PROCS_USED_TO_PRINT}"
      "${MAX_NUM_PROCESSORS_USED}"  "${TIMEOUT_USED}" )

    #
    # F.2) Write the cmake -P script
    #
  
    APPEND_STRING_VAR( TEST_SCRIPT_STR
      "\n"
      "SET(PROJECT_NAME ${PROJECT_NAME})\n"
      "\n"
      "SET(${PROJECT_NAME}_TRIBITS_DIR ${${PROJECT_NAME}_TRIBITS_DIR})\n"
      "\n"
      "SET(TEST_NAME ${TEST_NAME})\n"
      "\n"
      "SET(NUM_CMNDS ${NUM_CMNDS})\n"
      "\n"
      "SET(OVERALL_WORKING_DIRECTORY \"${PARSE_OVERALL_WORKING_DIRECTORY}\")\n"
      "\n"
      "SET(SKIP_CLEAN_OVERALL_WORKING_DIRECTORY \"${PARSE_SKIP_CLEAN_OVERALL_WORKING_DIRECTORY}\")\n"
      "\n"
      "SET(FAIL_FAST ${PARSE_FAIL_FAST})\n"
      "\n"
      "SET(SHOW_START_END_DATE_TIME ${${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME})\n"
      "\n"
      "SET(SHOW_MACHINE_LOAD ${${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST})\n"
      "\n"
      "SET(CATEGORIES ${PARSE_CATEGORIES})\n"
      "\n"
      "SET(PROCESSORS ${MAX_NUM_PROCESSORS_USED})\n"
      "\n"
      "SET(TIMEOUT ${TIMEOUT_USED})\n"
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
  
    IF (NOT TRIBITS_ADD_ADVANCED_TEST_SKIP_SCRIPT)
  
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("\nWriting file \"${TEST_SCRIPT_FILE}\" ...")
      ENDIF()
  
      FILE( WRITE "${TEST_SCRIPT_FILE}"
        "${TEST_SCRIPT_STR}" )
  
    ENDIF()


  ELSE()

    GLOBAL_SET(TRIBITS_ADD_ADVANCED_TEST_NUM_CMNDS "")

  ENDIF()

ENDFUNCTION()
