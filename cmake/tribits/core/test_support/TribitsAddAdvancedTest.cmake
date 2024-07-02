# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsCMakePolicies.cmake"  NO_POLICY_SCOPE)
include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsConstants.cmake")

set(tribitsAddAdvancedTestModuleDir "${CMAKE_CURRENT_LIST_DIR}")

include("${CMAKE_CURRENT_LIST_DIR}/TribitsAddAdvancedTestHelpers.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/../utils/TribitsPrintList.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../utils/AppendStringVar.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../utils/PrintVar.cmake")


# @FUNCTION: tribits_add_advanced_test()
#
# Function that creates an advanced test defined by stringing together one or
# more executable and/or command invocations that is run as a ``cmake -P``
# script with very flexible pass/fail criteria.
#
# Usage::
#
#   tribits_add_advanced_test(
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
#     [DISABLED <messageWhyDisabled>]
#     [FINAL_PASS_REGULAR_EXPRESSION "<regex>" |
#       FINAL_FAIL_REGULAR_EXPRESSION "<regex>"]
#     [ENVIRONMENT <var1>=<value1> <var2>=<value2> ...]
#     [TIMEOUT <maxSeconds>]
#     [LIST_SEPARATOR <sep>]
#     [ADDED_TEST_NAME_OUT <testName>]
#     )
#
# This function allows one to add a single CTest test that is actually a
# sequence of one or more separate commands strung together in some way to
# define the final pass/fail. One will want to use this function to add a test
# instead of `tribits_add_test()`_ when one needs to run more than one
# command, or one needs more sophisticated checking of the test result other
# than just grepping STDOUT (e.g. by running separate post-processing programs
# to examine output files).
#
# For more details on these arguments, see `TEST_<idx> EXEC/CMND Test Blocks
# and Arguments (tribits_add_advanced_test())`_.
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
#      [ARGS "<arg0>" "<arg1>" ... "<argn>"]
#      [MESSAGE "<message>"]
#      [WORKING_DIRECTORY <workingDir>]
#      [SKIP_CLEAN_WORKING_DIRECTORY]
#      [NUM_MPI_PROCS <numProcs>]
#      [NUM_TOTAL_CORES_USED <numTotalCoresUsed>]
#      [OUTPUT_FILE <outputFile>]
#      [NO_ECHO_OUTPUT]]
#      [PASS_ANY
#        | PASS_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...
#        | PASS_REGULAR_EXPRESSION_ALL "<regex0>" "<regex1>" ...
#        | STANDARD_PASS_OUTPUT ]
#      [FAIL_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...]
#      [ALWAYS_FAIL_ON_NONZERO_RETURN | ALWAYS_FAIL_ON_ZERO_RETURN]
#      [WILL_FAIL]
#
# For more information on these arguments, see `TEST_<idx> EXEC/CMND Test
# Blocks and Arguments (tribits_add_advanced_test())`_.
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
#   ``configure_file( ... COPYONLY)`` or build again when using
#   ``tribits_copy_files_to_binary_dir()`` in order to copy files).
#
# * When using ``OVERALL_WORKING_DIR TEST_NAME``, the test directory gets blow
#   away every time before it runs and therefore any old files are deleted
#   before the test gets run again (which avoids the problem of having a test
#   pass looking for the old files that will not be there when someone
#   configures and builds from scratch).
#
# For more information on these arguments, see `TEST_<idx>
# COPY_FILES_TO_TEST_DIR Test Blocks and Arguments
# (tribits_add_advanced_test())`_.
# 
# By default, each and every atomic ``TEST_<idx>`` block needs to pass (as
# defined in `Test case Pass/Fail (tribits_add_advanced_test())`_) in order
# for the overall test to pass.
#
# Finally, the test is only added if tests are enabled for the package
# (i.e. `${PACKAGE_NAME}_ENABLE_TESTS`_ ``= ON``) and if other criteria are
# met (see `Overall Arguments (tribits_add_advanced_test())`_).  (NOTE: A more
# efficient way to optionally enable tests is to put them in a ``test/``
# subdir and then include that subdir with `tribits_add_test_directories()`_.)
#
# *Sections:*
#
# * `Overall Arguments (tribits_add_advanced_test())`_
# * `TEST_<idx> EXEC/CMND Test Blocks and Arguments (tribits_add_advanced_test())`_
# * `TEST_<idx> COPY_FILES_TO_TEST_DIR Test Blocks and Arguments (tribits_add_advanced_test())`_
# * `Test case Pass/Fail (tribits_add_advanced_test())`_
# * `Overall Pass/Fail (tribits_add_advanced_test())`_
# * `Argument Parsing and Ordering (tribits_add_advanced_test())`_
# * `Implementation Details (tribits_add_advanced_test())`_
# * `Setting Additional Test Properties (tribits_add_advanced_test())`_
# * `Running multiple tests at the same time (tribits_add_advanced_test())`_
# * `Disabling Tests Externally (tribits_add_advanced_test())`_
# * `Debugging and Examining Test Generation (tribits_add_advanced_test())`_
# * `Using tribits_add_advanced_test() in non-TriBITS CMake projects`_
#
# .. _Overall Arguments (tribits_add_advanced_test()):
#
# **Overall Arguments (tribits_add_advanced_test())**
#
# Below, some of the overall arguments are described.  The rest of the overall
# arguments that control overall pass/fail are described in `Overall Pass/Fail
# (tribits_add_advanced_test())`_.  (NOTE: All of these arguments must be
# listed outside of the ``TEST_<idx>`` blocks, see `Argument Parsing and
# Ordering (tribits_add_advanced_test())`_).
#
#   ``<testNameBase>``
#
#     The base name of the test (which will have ``${PACKAGE_NAME}_``
#     prepended to the name, see <testName> below) that will be used to name
#     the output CMake script file as well as the CTest test name passed into
#     ``add_test()``.  This must be the first argument to this function.  The
#     name is allowed to contain '/' chars but these will be replaced with
#     '__' in the overall working directory name and the ctest -P script
#     (`Debugging and Examining Test Generation
#     (tribits_add_advanced_test())`_).
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
#     If specified, then no other tests will be allowed to run while this test
#     is running. See the ``RUN_SERIAL`` argument in the function
#     `tribits_add_test()`_ for more details.
#
#   ``COMM [serial] [mpi]``
#
#     If specified, selects if the test will be added in serial and/or MPI
#     mode.  See the ``COMM`` argument in the function `tribits_add_test()`_
#     for more details.
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
#     (tribits_add_advanced_test())`_).  **WARNING!** If just running a serial
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
#     See `tribits_add_test()`_ for more details.
#
#   ``HOST <host0> <host1> ...``
#
#     The list of hosts for which to enable the test (see
#     `tribits_add_test()`_).
#
#   ``XHOST <host0> <host1> ...``
#
#     The list of hosts for which **not** to enable the test (see
#     `tribits_add_test()`_).
#
#   ``HOSTTYPE <hosttype0> <hosttype1> ...``
#
#     The list of host types for which to enable the test (see
#     `tribits_add_test()`_).
#
#   ``XHOSTTYPE <hosttype0> <hosttype1> ...``
#
#     The list of host types for which **not** to enable the test (see
#     `tribits_add_test()`_).
#
#   ``EXCLUDE_IF_NOT_TRUE <varname0> <varname1> ...``
#
#     If specified, gives the names of CMake variables that must evaluate to
#     true, or the test will not be added (see `tribits_add_test()`_).
#
#   ``DISABLED <messageWhyDisabled>``
#
#     If ``<messageWhyDisabled>`` is non-empty and does not evaluate to FALSE
#     by CMake, then the test will be added by ctest but the ``DISABLED`` test
#     property will be set (see `tribits_add_test()`_).
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
#     to run before being timed-out and killed (see `Setting timeouts for
#     tests (tribits_add_test())`_).  This is for the full CTest test, not
#     individual ``TEST_<idx>`` commands!
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
#   ``ADDED_TEST_NAME_OUT <testName>``
#
#     If specified, then on output the variable ``<testName>`` will be set
#     with the name of the test passed to ``add_test()``.  Having this name
#     allows the calling ``CMakeLists.txt`` file access and set additional
#     test properties (see `Setting additional test properties
#     (tribits_add_advanced_test())`_).
#
# .. _TEST_<idx> EXEC/CMND Test Blocks and Arguments (tribits_add_advanced_test()):
#
#
# **TEST_<idx> EXEC/CMND Test Blocks and Arguments (tribits_add_advanced_test())**
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
#     the `tribits_add_test()`_ function (see `Determining the Executable or
#     Command to Run (tribits_add_test())`_).  If this is an MPI build, then
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
#     ``my_python_test.py`` with ``/usr/bin/env python`` at the top, you can't
#     just use::
#
#       CMND <path>/my_python_test.py ARGS "<arg0>" "<arg1>" ...
#
#     The same goes for Perl or any other scripting language.
#
#     Instead, you have to use::
#
#       CMND ${PYTHON_EXECUTABLE} ARGS <path>/my_python_test.py <arg0> <arg1> ...
#
#  ``ARGS "<arg0>" "<arg1>" ... "<argN>"``
#
#    The list of command-line arguments to pass to the ``CMND`` command or
#    ``EXEC`` executable.  Put each argument ``<argi>`` in quotes ``"<argi>"``
#    if it contains any spaces.  Also, of any of the individual arguments need
#    to contain semi-colons ``';'`` such as ``--my-arg=a;b a;c;d``, then pass
#    that quoted as ``"--my-arg=a<sep>b a<sep>c<sep>d"`` where ``<sep>``
#    matches the ``<sep>`` argument to the input ``LIST_SEPARATOR <sep>``.
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
#     given by relative ``<dir>``.  See `tribits_add_test()`_.
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
#     If that is not specified, then the value is taken from
#     ``${MPI_EXEC_DEFAULT_NUMPROCS}``.  For serial builds
#     (i.e. ``TPL_ENABLE_MPI=OFF``), passing in a value ``<numMpiProcs>`` >
#     ``1`` will cause the entire test to not be added.
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
#     be printed to STDOUT unless ``NO_ECHO_OUTPUT`` is passed as well.
#
#     NOTE: Contrary to CMake documentation for execute_process(), STDOUT and
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
# (see `Test case Pass/Fail (tribits_add_advanced_test())`_):
#
#   ``PASS_ANY``
#
#     If specified, the test command will be assumed to pass regardless of
#     the return value or any other output.  This would be used when a command
#     that is to follow will determine pass or fail based on output from this
#     command in some way.
#
#   ``PASS_REGULAR_EXPRESSION "<regex0>"  "<regex1>" ...``
#
#     If specified, the test command will be assumed to pass if it matches
#     **any** of the given regular expressions.  Otherwise, it is assumed to
#     fail.  TIPS: Replace ';' with '[;]' or CMake will interpret this as an
#     array element boundary.  To match '.', use '[.]'.
#
#   ``PASS_REGULAR_EXPRESSION_ALL "<regex0>" "<regex1>" ...``
#
#     If specified, the test command will be assumed to pass if the output
#     matches **all** of the provided regular expressions.  Note that this is not
#     a capability of raw ctest and represents an extension provided by
#     TriBITS.  NOTE: It is critical that you replace ';' with '[;]' or CMake
#     will interpret this as an array element boundary.
#
#   ``STANDARD_PASS_OUTPUT``
#
#     If specified, the test command will be assumed to pass if the string
#     expression "Final Result: PASSED" is found in the output for the test.
#     This as the result of directly passing in ``PASS_REGULAR_EXPRESSION
#     "End Result: TEST PASSED"``.
#
#   ``FAIL_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...``
#
#     If specified, the test command will be assumed to fail if it matches
#     **any** of the given regular expressions.  This will be applied and take
#     precedence over other above pass criteria.  For example, if even if
#     ``PASS_REGULAR_EXPRESSION`` or ``PASS_REGULAR_EXPRESSION_ALL`` match,
#     then the test will be marked as failed if any of the fail regexes match
#     the output.
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
# `Argument Parsing and Ordering (tribits_add_advanced_test())`_).
#
# **NOTE:** The current implementation limits the number of ``TEST_<idx>``
# blocks to just 20 (i.e. for ``<idx>=0...19``).  If more test blocks are
# added (e.g. ``TEST_20``), then an fatal error message will be printed and
# processing will end.  To increase this max in a local scope, call::
#
#   set(TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS <larger-num>)
#
# where ``<larger-num> > 20``.  This can be set in any scope in any
# ``CMakeLists.txt`` file or inside of a function and it will impact all of
# the future calls to ``tribits_add_advanced_test()`` in that scope.
#
# .. _TEST_<idx> COPY_FILES_TO_TEST_DIR Test Blocks and Arguments (tribits_add_advanced_test()):
#
# **TEST_<idx> COPY_FILES_TO_TEST_DIR Test Blocks and Arguments (tribits_add_advanced_test())**
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
# .. _Test case Pass/Fail (tribits_add_advanced_test()):
#
# **Test case Pass/Fail (tribits_add_advanced_test())**
#
# The logic given below can be used to determine pass/fail criteria for a test
# case both based on what is printed in the test output **and** the return
# code for the test block command.  Raw CTest, as of version 3.23, does not
# allow that.  With raw CTest, one can only set pass/fail criteria based the
# test output **or** the return code, but not both.  This make
# `tribits_add_advanced_test()`_ more attractive to use than
# `tribits_add_test()`_ or raw ``add_test()`` in cases where it is important
# to check both.
#
# The logic for how pass/fail for a ``TEST_<IDX>`` ``EXEC`` or ``CMND`` case
# is applied is given by::
#
#   # A) Apply first set of pass/fail logic
#   TEST_CASE_PASSED = FALSE
#   If PASS_ANY specified:
#     TEST_CASE_PASSED = TRUE
#   Else If PASS_REGULAR_EXPRESSION is specified:
#     For each "<regexi>" in PASS_REGULAR_EXPRESSION:
#       If "<regexi>" matches STDOUT:
#         TEST_CASE_PASSED = TRUE
#       Endif
#     Endforeach
#   Else if PASS_REGULAR_EXPRESSION_ALL specified:
#     TEST_CASE_PASSED = TRUE
#     For each "<regexi>" in PASS_REGULAR_EXPRESSION_ALL:
#       If "<regexi>" does not match STDOUT:
#         TEST_CASE_PASSED = FALSE
#       Endif
#     Endforeach
#   Else
#     If command return code == 0:
#       TEST_CASE_PASSED = TRUE
#     Endif
#   Endif
#
#   # B) Check for failing regex matching?
#   If FAIL_REGULAR_EXPRESSION specified:
#     For each "<regexi>" in FAIL_REGULAR_EXPRESSION:
#       If "<regexi>" matches STDOUT:
#         TEST_CASE_PASSED = FALSE
#       Endif
#     Endforeach
#   Endif
#
#   # C) Check for return code always 0 or !=0?
#   If ALWAYS_FAIL_ON_NONZERO_RETURN specified and return code != 0:
#     TEST_CASE_PASSED = FALSE
#   Else If ALWAYS_FAIL_ON_ZERO_RETURN specified and return code == 0:
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
# Note that the above is the exact same logic that CTest uses to determine
# pass/fail w.r.t. to the CTest properties ``PASS_REGULAR_EXPRESSION``,
# ``FAIL_REGULAR_EXPRESSION`` and ``WILL_FAIL``.  (It is just that raw
# CMake/CTest, as of version 3.23, does not support any pass/fail criteria
# like ``PASS_REGULAR_EXPRESSION_ALL`` or
# ``ALWAYS_FAIL_ON_NONZERO_RETURN``/``ALWAYS_FAIL_ON_ZERO_RETURN``.)
#
# .. _Overall Pass/Fail (tribits_add_advanced_test()):
#
# **Overall Pass/Fail (tribits_add_advanced_test())**
#
# By default, the overall test will be assumed to pass if it prints::
#
#   "OVERALL FINAL RESULT: TEST PASSED (<testName>)"
#
# However, this can be changed by setting one of the following optional arguments:
#
#   ``FINAL_PASS_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...``
#
#     If specified, the test will be assumed to pass if the output matches
#     **any** of the provided regular expressions ``<regexi>``.  (Sets the
#     CTest property ``PASS_REGULAR_EXPRESSION`` for the overall test.)
#
#   ``FINAL_FAIL_REGULAR_EXPRESSION "<regex0>" "<regex1>" ...``
#
#     If specified, the test will be assumed to fail if the output matches
#     **any** of the provided regular expressions ``<regexi>`` regardless if
#     other criteria would have the test passing.  (Sets the CTest property
#     ``FAIL_REGULAR_EXPRESSION`` for the overall test.)
#
# **NOTE:** It is **not** recommended to set ``FINAL_PASS_REGULAR_EXPRESSION``
# or ``FINAL_FAIL_REGULAR_EXPRESSION`` directly, but instead to determine
# pass/fail for each test case individually as described in `TEST_<idx>
# EXEC/CMND Test Blocks and Arguments (tribits_add_advanced_test())`_ and
# `Test case Pass/Fail (tribits_add_advanced_test())`_.  Otherwise, the test
# will confuse most people and the output behavior will seem very strange.
#
# .. _Argument Parsing and Ordering (tribits_add_advanced_test()):
#
# **Argument Parsing and Ordering (tribits_add_advanced_test())**
#
# The basic tool used for parsing the arguments to this function is the
# command ``cmake_parse_arguments()`` which has a certain set of behaviors.
# The parsing using ``cmake_parse_arguments()`` is actually done in two
# phases.  There is a top-level parsing of the "overall" arguments listed in
# `Overall Arguments (tribits_add_advanced_test())`_ that also pulls out the
# test blocks.  Then there is a second level of parsing using
# ``cmake_parse_arguments()`` for each of the ``TEST_<idx>`` blocks.  Because
# of this usage, there are a few restrictions that one needs to be aware of
# when using ``tribits_add_advanced_test()``.  This short sections tries to
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
# .. _Implementation Details (tribits_add_advanced_test()):
#
# **Implementation Details (tribits_add_advanced_test())**
#
# Since raw CTest does not support the features provided by this function, the
# way an advanced test is implemented is that a ``cmake -P`` script with the
# name ``<testName>.cmake`` (with any '/' replaced with '__') gets created in
# the current binary directory that then gets added to CTest using::
#
#   add_test(<testName> cmake [other options] -P <testName>.cmake)
#
# This ``cmake -P`` script then runs the various test cases and checks the
# pass/fail for each case to determine overall pass/fail and implement other
# functionality described above.
#
# .. _Setting Additional Test Properties (tribits_add_advanced_test()):
#
# **Setting Additional Test Properties (tribits_add_advanced_test())**
#
# After this function returns, if the test gets added using ``add_test()``,
# then additional properties can be set and changed using
# ``set_tests_properties(<testName> ...)``, where ``<testName>`` is returned
# using the ``ADDED_TEST_NAME_OUT <testName>`` argument.  Therefore, any tests
# properties that are not directly supported by this function and passed
# through the argument list to this wrapper function can be set in the outer
# ``CMakeLists.txt`` file after the call to ``tribits_add_advanced_test()``.
# For example::
#
#   tribits_add_advanced_test_test( someTest ...
#     ADDED_TEST_NAME_OUT  someTest_TEST_NAME )
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
# The following built-in CTest test properties are set through `Overall
# Arguments (tribits_add_advanced_test())`_ or are otherwise automatically set
# by this function and should **NOT** be overridden by direct calls to
# ``set_tests_properties()``: ``ENVIRONMENT``, ``FAIL_REGULAR_EXPRESSION``,
# ``LABELS``, ``PASS_REGULAR_EXPRESSION``, ``RUN_SERIAL``, ``TIMEOUT``,
# ``WILL_FAIL``, and ``WORKING_DIRECTORY``.
#
# However, generally, other built-in CTest test properties can be set after
# the test is added like show above.  Examples of test properties that can be
# set using direct calls to ``set_tests_properties()`` include
# ``ATTACHED_FILES``, ``ATTACHED_FILES_ON_FAIL``, ``COST``, ``DEPENDS``,
# ``MEASUREMENT``, and ``RESOURCE_LOCK``.
#
# For example, one can set a dependency between two tests using::
#
#   tribits_add_advanced_test_test( test_a [...]
#      ADDED_TEST_NAME_OUT  test_a_TEST_NAME )
#   
#   tribits_add_advanced_test_test( test_b [...]
#      ADDED_TEST_NAME_OUT  test_z_TEST_NAME )
#   
#   if (test_a_TEST_NAME AND test_b_TEST_NAME)
#     set_tests_properties(${test_b_TEST_NAME}
#       PROPERTIES DEPENDS ${test_a_TEST_NAME})
#   endif()
#
# This ensures that test ``test_b`` will always be run after ``test_a`` if
# both tests are run by CTest.
#
# .. _Running multiple tests at the same time (tribits_add_advanced_test()):
#
# **Running multiple tests at the same time (tribits_add_advanced_test())**
#
# Just as with `tribits_add_test()`_, setting ``NUM_MPI_PROCS <numProcs>`` or
# ``OVERALL_NUM_MPI_PROCS <numOverallProcs>`` or ``NUM_TOTAL_CORES_USED
# <numTotalCoresUsed>`` or ``OVERALL_NUM_TOTAL_CORES_USED
# <overallNumTotalCoresUsed>`` will set the ``PROCESSORS`` CTest property to
# allow CTest to schedule and run multiple tests at the same time when ``'ctest
# -j<N>'`` is used (see `Running multiple tests at the same time
# (tribits_add_test())`_).
#
# .. _Disabling Tests Externally (tribits_add_advanced_test()):
#
# **Disabling Tests Externally (tribits_add_advanced_test())**
#
# The test can be disabled externally by setting the CMake cache variable
# ``<testName>_DISABLE=TRUE``.  This allows tests to be disabled on a
# case-by-case basis.  The name ``<testName>`` must be the *exact* name
# that shows up in ``ctest -N`` when running the test.
#
# .. _Debugging and Examining Test Generation (tribits_add_advanced_test()):
#
# **Debugging and Examining Test Generation (tribits_add_advanced_test())**
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
# directory (see `Implementation Details (tribits_add_advanced_test())`_) and
# the generated ``CTestTestfile.cmake`` file that should list this test case.
#
# .. _Using tribits_add_advanced_test() in non-TriBITS CMake projects:
#
# **Using tribits_add_advanced_test() in non-TriBITS CMake projects**
#
# The function ``tribits_add_advanced_test()`` can be used to add tests in
# non-TriBITS projects.  To do so, one just needs to set the variables
# ``${PROJECT_NAME}_ENABLE_TESTS=TRUE`` and ``${PROJECT_NAME}_TRIBITS_DIR``
# (pointing to the TriBITS location).  For example, a valid project can be a
# simple as::
#
#   cmake_minimum_required(VERSION 3.23.0)
#   set(PROJECT_NAME TAATDriver)
#   project(${PROJECT_NAME} NONE)
#   set(${PROJECT_NAME}_TRACE_ADD_TEST TRUE)
#   set(${PROJECT_NAME}_TRIBITS_DIR ""  CACHE FILEPATH
#     "Location of TriBITS to use." ) 
#   set(PACKAGE_NAME ${PROJECT_NAME})
#   set(${PACKAGE_NAME}_ENABLE_TESTS TRUE)
#   include("${${PROJECT_NAME}_TRIBITS_DIR}/core/test_support/TribitsAddAdvancedTest.cmake")
#   include(CTest)
#   enable_testing()
#   
#   tribits_add_advanced_test( HelloWorld
#     OVERALL_WORKING_DIRECTORY TEST_NAME
#     TEST_0 CMND echo ARGS "Hello World!"
#       PASS_REGULAR_EXPRESIOIN "Hello World"
#     )
#
# Above, one can replace::
#
#   include("${${PROJECT_NAME}_TRIBITS_DIR}/core/test_support/TribitsAddAdvancedTest.cmake")
#
# with::
#
#   list(PREPEND CMAKE_MODULE_PATH "${${PROJECT_NAME}_TRIBITS_DIR}/core/test_support")
#   include(TribitsAddAdvancedTest)
#
# and it will have the same effect.
#
function(tribits_add_advanced_test TEST_NAME_IN)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nPACKAGE_ADD_ADVANCED_TEST: ${TEST_NAME_IN}\n")
  endif()

  tribits_set_tribits_package_name()

  global_set(TRIBITS_SET_TEST_PROPERTIES_INPUT)
  global_set(MESSAGE_WRAPPER_INPUT)

  # Set the full TEST_NAME
  set(TEST_NAME ${PACKAGE_NAME}_${TEST_NAME_IN})

  #
  # A) Parse the overall arguments and figure out how many tests
  # commands we will have
  #

  # Set maximum number of TEST_<idx> blocks
  tribits_add_advanced_test_max_num_test_cmnd_idx_compute()
  set(MAX_NUM_TEST_CMND_IDX ${TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX})

  set(TEST_IDX_LIST "")
  foreach( TEST_CMND_IDX RANGE ${MAX_NUM_TEST_CMND_IDX})
    list( APPEND TEST_IDX_LIST TEST_${TEST_CMND_IDX} )
  endforeach()

  set(optionsList  FAIL_FAST  RUN_SERIAL  SKIP_CLEAN_OVERALL_WORKING_DIRECTORY)

  set(oneValueKeywordsList  DISABLED)

  set(multiValueKeywordsList
    ${TEST_IDX_LIST}  OVERALL_WORKING_DIRECTORY
    LIST_SEPARATOR
    OVERALL_NUM_MPI_PROCS  OVERALL_NUM_TOTAL_CORES_USED
    CATEGORIES  COMM  HOST  XHOST  HOSTTYPE  XHOSTTYPE  EXCLUDE_IF_NOT_TRUE
    FINAL_PASS_REGULAR_EXPRESSION  FINAL_FAIL_REGULAR_EXPRESSION
    TIMEOUT  ENVIRONMENT  KEYWORDS
    ADDED_TEST_NAME_OUT
    )

  cmake_parse_arguments(
    PARSE_ARGV 1  # NOTE: One named argument to skip over
    PARSE  # prefix
    "${optionsList}"
    "${oneValueKeywordsList}"
    "${multiValueKeywordsList}"
    )

  tribits_check_for_unparsed_arguments()
  tribits_add_advanced_test_check_exceed_max_num_test_blocks()

  if(PARSE_ADDED_TEST_NAME_OUT)
    set(${PARSE_ADDED_TEST_NAME_OUT} "" PARENT_SCOPE )
  endif()

  # Set the name of the cmake -P script.
  string(REPLACE "/" "__" NORMALIZED_TEST_NAME "${TEST_NAME}")
  set(TEST_SCRIPT_FILE_NAME "${NORMALIZED_TEST_NAME}.cmake")

  # Set the relative overall working directory and abs working directory
  if (PARSE_OVERALL_WORKING_DIRECTORY)
    if ("${PARSE_OVERALL_WORKING_DIRECTORY}" STREQUAL "TEST_NAME")
      set(PARSE_OVERALL_WORKING_DIRECTORY ${NORMALIZED_TEST_NAME})
    endif()
   # Test will run in created working subdir
   set(ABS_OVERALL_WORKING_DIRECTORY
     ${CMAKE_CURRENT_BINARY_DIR}/${PARSE_OVERALL_WORKING_DIRECTORY})
  else()
    # Test runs in current binary directory (not a good idea!) 
    set(ABS_OVERALL_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()

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

  tribits_add_test_query_disable(DISABLE_THIS_TEST ${TEST_NAME})
  if (DISABLE_THIS_TEST)
    return()
  endif()

  tribits_set_run_serial(${TEST_NAME} "${PARSE_RUN_SERIAL}"
    SET_RUN_SERIAL)

  tribits_set_disabled_and_msg(${TEST_NAME} "${PARSE_DISABLED}"
    SET_DISABLED_AND_MSG)  # Adds the test but sets DISABLED test prop!

  #
  # C) Determine if we will add the serial or MPI tests based on input COMM
  # and TPL_ENABLE_MPI
  #

  tribits_process_comm_args(ADD_SERIAL_TEST  ADD_MPI_TEST  ${PARSE_COMM})
  if (NOT ADD_SERIAL_TEST AND NOT ADD_MPI_TEST)
    return()
  endif()

  #
  # D) Build the test script
  #

  set(ADD_THE_TEST TRUE)

  set(TEST_SCRIPT_STR "")

  string(APPEND  TEST_SCRIPT_STR
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
    "set( TEST_NAME ${TEST_NAME} )\n"
    "set( LIST_SEPARATOR \"${PARSE_LIST_SEPARATOR}\" )\n"
    )

  # Loop through each test case

  set(NUM_CMNDS 0)
  set(TEST_EXE_LIST "")

  if (PARSE_OVERALL_NUM_MPI_PROCS  AND  PARSE_OVERALL_NUM_TOTAL_CORES_USED)
    if (PARSE_OVERALL_NUM_MPI_PROCS  GREATER  PARSE_OVERALL_NUM_TOTAL_CORES_USED)
      message_wrapper(FATAL_ERROR
        "ERROR: ${TEST_NAME}: OVERALL_NUM_MPI_PROCS='${PARSE_OVERALL_NUM_MPI_PROCS}' > OVERALL_NUM_TOTAL_CORES_USED='${PARSE_OVERALL_NUM_TOTAL_CORES_USED}' not allowed!")
      return()
    endif()
  endif()

  # ToDo: Assert that OVERALL_NUM_TOTAL_CORES_USED >= OVERALL_NUM_MPI_PROCS

  if (PARSE_OVERALL_NUM_MPI_PROCS)
    set(MAX_NUM_MPI_PROCS_USED ${PARSE_OVERALL_NUM_MPI_PROCS})
    set(MAX_NUM_PROCESSORS_USED ${PARSE_OVERALL_NUM_MPI_PROCS})
  else()
    set(MAX_NUM_MPI_PROCS_USED ${PARSE_OVERALL_NUM_MPI_PROCS})
    set(MAX_NUM_PROCESSORS_USED 1)
  endif()

  set(HAS_AT_LEAST_ONE_EXEC FALSE)

  foreach( TEST_CMND_IDX RANGE ${MAX_NUM_TEST_CMND_IDX} )

    if (NOT PARSE_TEST_${TEST_CMND_IDX} )
      break()
    endif()

    if (NOT  ADD_THE_TEST)
      break()
    endif()

    math( EXPR NUM_CMNDS ${NUM_CMNDS}+1 )

    # Parse the test command case

    # Search to see if we are copying files or not for this TEST_<IDX> block ...

    set(PARSE_COPY_FILES_TO_TEST_DIR)
    set(COPY_FILES_TO_TEST_DIR_IDX FALSE)
    foreach(PARSE_TEST_IDX_ARGS ${PARSE_TEST_${TEST_CMND_IDX}})
      if (PARSE_TEST_IDX_ARGS STREQUAL "COPY_FILES_TO_TEST_DIR")
        set(COPY_FILES_TO_TEST_DIR_IDX TRUE)
      endif()
    endforeach()

    if (COPY_FILES_TO_TEST_DIR_IDX)

      # Do a special parse just for TEST_<IDX> blocks of type
      # COPY_FILES_TO_TEST_DIR

      cmake_parse_arguments(
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
      tribits_check_for_unparsed_arguments()
      tribits_assert_parse_arg_one_or_more_values(PARSE COPY_FILES_TO_TEST_DIR)
      tribits_assert_parse_arg_zero_or_one_value(PARSE SOURCE_DIR)
      tribits_assert_parse_arg_zero_or_one_value(PARSE DEST_DIR)

      set(PARSE_EXEC "")
      set(PARSE_CMND "")

    else()

      # Parse TEST_<IDX> block args for types EXEC and CMND

      set(testBlockOptionsList  NOEXEPREFIX  NOEXESUFFIX  NO_ECHO_OUTPUT  PASS_ANY
        STANDARD_PASS_OUTPUT  ALWAYS_FAIL_ON_NONZERO_RETURN  ALWAYS_FAIL_ON_ZERO_RETURN
        WILL_FAIL  ADD_DIR_TO_NAME  SKIP_CLEAN_WORKING_DIRECTORY
        )

      set(testBlockMultiValueKeywordsList  EXEC  CMND  ARGS  DIRECTORY  MESSAGE
        WORKING_DIRECTORY  OUTPUT_FILE  NUM_MPI_PROCS  NUM_TOTAL_CORES_USED
        PASS_REGULAR_EXPRESSION_ALL  FAIL_REGULAR_EXPRESSION  PASS_REGULAR_EXPRESSION
        )

      cmake_parse_arguments(
         PARSE  #prefix
         "${testBlockOptionsList}"
         ""     # one_value_keywords
         "${testBlockMultiValueKeywordsList}"
         ${PARSE_TEST_${TEST_CMND_IDX}}
         )

      tribits_check_for_unparsed_arguments(PARSE) # ToDo: Use a different prefix!

    endif()

    #
    # Set up the command that will be written into the cmake -P *.cmake file
    #

    set(ARGS_STR "${PARSE_ARGS}")

    if (PARSE_EXEC)

      #
      # This is an EXEC test block
      #

      set(HAS_AT_LEAST_ONE_EXEC TRUE)

      list( LENGTH PARSE_EXEC PARSE_EXEC_LEN )
      if (NOT PARSE_EXEC_LEN EQUAL 1)
        message(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} EXEC = '${PARSE_EXEC}'"
          " must be a single name.  To add arguments use ARGS <arg1> <arg2> ...." )
      endif()

      tribits_add_test_get_exe_binary_name( "${PARSE_EXEC}"
        ${PARSE_NOEXEPREFIX} ${PARSE_NOEXESUFFIX}
        ${PARSE_ADD_DIR_TO_NAME} EXE_BINARY_NAME )

      tribits_add_test_adjust_directory( ${EXE_BINARY_NAME} "${PARSE_DIRECTORY}"
        EXECUTABLE_PATH)

      if (PARSE_NUM_MPI_PROCS)
        set(NUM_MPI_PROC_VAR_NAME "NUM_MPI_PROCS")
      else()
        set(PARSE_NUM_MPI_PROCS ${PARSE_OVERALL_NUM_MPI_PROCS})
        set(NUM_MPI_PROC_VAR_NAME "OVERALL_NUM_MPI_PROCS")
      endif()

      tribits_add_test_get_num_procs_used("${PARSE_NUM_MPI_PROCS}"
        "${NUM_MPI_PROC_VAR_NAME}"  NUM_PROCS_USED  NUM_PROCS_USED_NAME)

      if (NUM_PROCS_USED LESS 0)
        set(ADD_THE_TEST FALSE)
      elseif (NUM_PROCS_USED  GREATER  MAX_NUM_MPI_PROCS_USED)
        set(MAX_NUM_MPI_PROCS_USED  ${NUM_PROCS_USED})
      endif()

      if (PARSE_NUM_TOTAL_CORES_USED)
        set(NUM_TOTAL_CORES_USED  ${PARSE_NUM_TOTAL_CORES_USED})
        set(NUM_TOTAL_CORES_USED_NAME  "NUM_TOTAL_CORES_USED")
      else()
        set(NUM_TOTAL_CORES_USED  ${PARSE_OVERALL_NUM_TOTAL_CORES_USED})
        set(NUM_TOTAL_CORES_USED_NAME  "OVERALL_NUM_TOTAL_CORES_USED")
      endif()

      tribits_add_test_get_num_total_cores_used("${TEST_NAME}${MPI_NAME_POSTFIX}"
        "${NUM_TOTAL_CORES_USED}"  "${NUM_TOTAL_CORES_USED_NAME}"
        "${NUM_PROCS_USED}"  "${NUM_PROCS_USED_NAME}"
        NUM_TOTAL_CORES_USED  SKIP_TEST)
      if (SKIP_TEST)
        set(ADD_THE_TEST FALSE)
      endif()

      if (NUM_TOTAL_CORES_USED  GREATER  MAX_NUM_PROCESSORS_USED)
        set(MAX_NUM_PROCESSORS_USED  ${NUM_TOTAL_CORES_USED})
      endif()

      if(ADD_THE_TEST)
        list(APPEND TEST_EXE_LIST ${EXECUTABLE_PATH})
      endif()

      tribits_add_test_get_test_cmnd_array( TEST_CMND_ARRAY
        "${EXECUTABLE_PATH}" "${NUM_PROCS_USED}" ${ARGS_STR} )

    elseif (PARSE_CMND)

      #
      # This is a COMMAND test block
      #

      list( LENGTH PARSE_CMND PARSE_CMND_LEN )
      if (NOT PARSE_CMND_LEN EQUAL 1)
        message_wrapper(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} CMND = '${PARSE_CMND}'"
          " must be a single command.  To add arguments use ARGS <arg1> <arg2> ...." )
      endif()

      if (PARSE_NUM_TOTAL_CORES_USED)
        set(NUM_TOTAL_CORES_USED  ${PARSE_NUM_TOTAL_CORES_USED})
        set(NUM_TOTAL_CORES_USED_NAME  "NUM_TOTAL_CORES_USED")
      else()
        set(NUM_TOTAL_CORES_USED  ${PARSE_OVERALL_NUM_TOTAL_CORES_USED})
        set(NUM_TOTAL_CORES_USED_NAME  "OVERALL_NUM_TOTAL_CORES_USED")
      endif()

      tribits_add_test_get_num_total_cores_used("${TEST_NAME}${MPI_NAME_POSTFIX}"
        "${NUM_TOTAL_CORES_USED}"  "${NUM_TOTAL_CORES_USED_NAME}"
        "1"  "DUMMY_NUM_MPI_PROCS"  # Never be printed
        NUM_TOTAL_CORES_USED  SKIP_TEST)
      if (SKIP_TEST)
        set(ADD_THE_TEST FALSE)
      endif()

      if (NUM_TOTAL_CORES_USED  GREATER  MAX_NUM_PROCESSORS_USED)
        set(MAX_NUM_PROCESSORS_USED  ${NUM_TOTAL_CORES_USED})
      endif()

      if(ADD_THE_TEST)
        if (NOT TRIBITS_ADD_TEST_ADD_TEST_UNITTEST)
          find_program(CMND_PATH ${PARSE_CMND})
        else()
          set(CMND_PATH ${PARSE_CMND})
        endif()
        list(APPEND TEST_EXE_LIST ${CMND_PATH})
      endif()

      set( TEST_CMND_ARRAY ${PARSE_CMND} ${ARGS_STR} )

    elseif (PARSE_COPY_FILES_TO_TEST_DIR)

      #
      # This is a COPY_FLES_TO_TEST_DIR block
      #

      # FILES_TO_COPY_COMMA_SEP
      set(FILES_TO_COPY_COMMA_SEP "${PARSE_COPY_FILES_TO_TEST_DIR}")
      string(REPLACE ";" "," FILES_TO_COPY_COMMA_SEP
        "${FILES_TO_COPY_COMMA_SEP}" )
      # NOTE: Above, we have to replace ';' with ',' or the lower commands
      # string(APPEND ) will replace ';' with ''.  This is *not* what we
      # want.  In DriveAdvancedTest.cmake, we will replace the ',' with ';'
      # again :-)  

      # SOURCE_DIR
      if (PARSE_SOURCE_DIR)
        if (IS_ABSOLUTE "${PARSE_SOURCE_DIR}")
          set(COPY_FILES_TO_TEST_DIR_SOURCE_DIR
            "${PARSE_SOURCE_DIR}")
        else()
          set(COPY_FILES_TO_TEST_DIR_SOURCE_DIR
            "${CMAKE_CURRENT_SOURCE_DIR}/${PARSE_SOURCE_DIR}")
        endif()
      else()
        set(COPY_FILES_TO_TEST_DIR_SOURCE_DIR
          "${CMAKE_CURRENT_SOURCE_DIR}")
      endif()

      # DEST_DIR
      if (PARSE_DEST_DIR)
        if (IS_ABSOLUTE "${PARSE_DEST_DIR}")
          set(COPY_FILES_TO_TEST_DIR_DEST_DIR
            "${PARSE_DEST_DIR}")
        else()
          set(COPY_FILES_TO_TEST_DIR_DEST_DIR
            "${ABS_OVERALL_WORKING_DIRECTORY}/${PARSE_DEST_DIR}")
        endif()
      else()
        set(COPY_FILES_TO_TEST_DIR_DEST_DIR
          "${ABS_OVERALL_WORKING_DIRECTORY}")
      endif()

    else()

      message( FATAL_ERROR
        "Must have EXEC, CMND, or COPY_FILES_TO_TEST_DIR for TEST_${TEST_CMND_IDX}" )

    endif()

    #
    # Write parts for this TEST_<IDX> block to TEST_SCRIPT_STR
    #

    if (PARSE_COPY_FILES_TO_TEST_DIR)

      # Write the vars for COPY_FILES_TO_TEST_DIR 
  
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_COPY_FILES_TO_TEST_DIR"
        " \"${FILES_TO_COPY_COMMA_SEP}\")\n"
        )
      if (TRIBITS_ADD_ADVANCED_TEST_UNITTEST)
        global_set(TRIBITS_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX}
          "${TEST_CMND_STR}" )
      endif()
  
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_SOURCE_DIR"
        " \"${COPY_FILES_TO_TEST_DIR_SOURCE_DIR}\")\n"
        )
  
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_DEST_DIR"
        " \"${COPY_FILES_TO_TEST_DIR_DEST_DIR}\")\n"
        )

    else()

      # Write the command to be run for EXEC and CMND blocks ...

      tribits_join_exec_process_set_args( TEST_CMND_STR "${TEST_CMND_ARRAY}" )

      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_CMND ${TEST_CMND_STR} )\n"
        )
      if (TRIBITS_ADD_ADVANCED_TEST_UNITTEST)
        global_set(TRIBITS_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX}
          "${TEST_CMND_STR}" )
      endif()

    endif()

    if (PARSE_MESSAGE)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_MESSAGE \"${PARSE_MESSAGE}\" )\n"
        )
    endif()

    if (PARSE_WORKING_DIRECTORY)
      if ("${PARSE_WORKING_DIRECTORY}" STREQUAL "TEST_NAME")
        set(PARSE_WORKING_DIRECTORY ${TEST_NAME})
      endif()
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_WORKING_DIRECTORY \"${PARSE_WORKING_DIRECTORY}\" )\n"
         )
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_SKIP_CLEAN_WORKING_DIRECTORY ${PARSE_SKIP_CLEAN_WORKING_DIRECTORY} )\n"
        )
    endif()

    if (PARSE_OUTPUT_FILE)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_OUTPUT_FILE \"${PARSE_OUTPUT_FILE}\" )\n"
        )
    endif()

    if (PARSE_NO_ECHO_OUTPUT)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_NO_ECHO_OUTPUT \"${PARSE_NO_ECHO_OUTPUT}\" )\n"
        )
    endif()

    # Set up pass/fail

    if (PARSE_PASS_ANY)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_PASS_ANY TRUE )\n"
        )
    elseif (PARSE_STANDARD_PASS_OUTPUT)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_PASS_REGULAR_EXPRESSION \"End Result: TEST PASSED\" )\n"
        )
    elseif (PARSE_PASS_REGULAR_EXPRESSION)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_PASS_REGULAR_EXPRESSION \"${PARSE_PASS_REGULAR_EXPRESSION}\" )\n"
        )
    elseif (PARSE_PASS_REGULAR_EXPRESSION_ALL)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_PASS_REGULAR_EXPRESSION_ALL "
        )
      foreach(REGEX_STR ${PARSE_PASS_REGULAR_EXPRESSION_ALL})
        string(APPEND  TEST_SCRIPT_STR
          "\"${REGEX_STR}\" "
          )
      endforeach()
      string(APPEND  TEST_SCRIPT_STR
        ")\n"
        )
    endif()

    if (PARSE_FAIL_REGULAR_EXPRESSION)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_FAIL_REGULAR_EXPRESSION \"${PARSE_FAIL_REGULAR_EXPRESSION}\" )\n"
        )
    endif()

    if (PARSE_ALWAYS_FAIL_ON_NONZERO_RETURN)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_ALWAYS_FAIL_ON_NONZERO_RETURN TRUE )\n"
        )
    endif()

    if (PARSE_ALWAYS_FAIL_ON_ZERO_RETURN)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_ALWAYS_FAIL_ON_ZERO_RETURN TRUE )\n"
        )
    endif()

    if (PARSE_WILL_FAIL)
      string(APPEND  TEST_SCRIPT_STR
        "\n"
        "set( TEST_${TEST_CMND_IDX}_WILL_FAIL TRUE )\n"
        )
    endif()

  endforeach()

  # ToDo: Verify that TEST_${MAX_NUM_TEST_CMND_IDX}+1 does *not* exist!

  #
  # F) Set the CTest test to run the new script
  #

  if (ADD_THE_TEST)

    #
    # F.1) Call add_test() and set the test properties
    #
  
    set(TEST_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_SCRIPT_FILE_NAME}")

    if(NOT TRIBITS_ADD_TEST_ADD_TEST_UNITTEST)
      # Tell CTest to run our script for this test.  Pass the test-type
      # configuration name to the script in the TEST_CONFIG variable.
      add_test( NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND} "-DTEST_CONFIG=\${CTEST_CONFIGURATION_TYPE}"
        -P "${TEST_SCRIPT_FILE}")
    endif()

    if(PARSE_ADDED_TEST_NAME_OUT)
      set(${PARSE_ADDED_TEST_NAME_OUT} ${TEST_NAME} PARENT_SCOPE )
    endif()

    list(REMOVE_DUPLICATES TEST_EXE_LIST)
    tribits_set_test_property(${TEST_NAME} PROPERTY REQUIRED_FILES ${TEST_EXE_LIST})

    if(SET_RUN_SERIAL)
      tribits_set_tests_properties(${TEST_NAME} PROPERTIES RUN_SERIAL ON)
    endif()

    tribits_private_add_test_add_label_and_keywords(${TEST_NAME})

    #This if clause will set the number of PROCESSORS to reserve during testing
    #to the number requested for the test.
    if(MAX_NUM_PROCESSORS_USED)
      tribits_set_tests_properties(${TEST_NAME} PROPERTIES
        PROCESSORS "${MAX_NUM_PROCESSORS_USED}")
    endif()

    tribits_private_add_test_add_environment_and_resource(${TEST_NAME}
      ${MAX_NUM_PROCESSORS_USED})

    if (SET_DISABLED_AND_MSG)
      tribits_set_tests_properties(${TEST_NAME} PROPERTIES DISABLED ON)
    endif()

    if (PARSE_FINAL_PASS_REGULAR_EXPRESSION)
      tribits_set_tests_properties( ${TEST_NAME} PROPERTIES
        PASS_REGULAR_EXPRESSION "${PARSE_FINAL_PASS_REGULAR_EXPRESSION}" )
    elseif (PARSE_FINAL_FAIL_REGULAR_EXPRESSION)
      tribits_set_tests_properties( ${TEST_NAME} PROPERTIES
        FAIL_REGULAR_EXPRESSION "${PARSE_FINAL_FAIL_REGULAR_EXPRESSION}" )
    else()
      tribits_set_tests_properties( ${TEST_NAME} PROPERTIES
        PASS_REGULAR_EXPRESSION
        "OVERALL FINAL RESULT: TEST PASSED .${TEST_NAME}." )
    endif()

    tribits_private_add_test_set_timeout(${TEST_NAME}  TIMEOUT_USED)

    tribits_private_add_test_set_environment(${TEST_NAME})

    if (TPL_ENABLE_MPI AND HAS_AT_LEAST_ONE_EXEC)
      set(MAX_NUM_MPI_PROCS_USED_TO_PRINT  ${MAX_NUM_MPI_PROCS_USED})
    else()
      set(MAX_NUM_MPI_PROCS_USED_TO_PRINT "")
    endif()

    tribits_private_add_test_print_added(${TEST_NAME}
      "${PARSE_CATEGORIES}"  "${MAX_NUM_MPI_PROCS_USED_TO_PRINT}"
      "${MAX_NUM_PROCESSORS_USED}"  "${TIMEOUT_USED}"
      "${SET_RUN_SERIAL}" "${SET_DISABLED_AND_MSG}")

    #
    # F.2) Write the cmake -P script
    #
  
    set(coreUtilsDir "${tribitsAddAdvancedTestModuleDir}/../utils")
    cmake_path(NORMAL_PATH coreUtilsDir) 
    string(APPEND  TEST_SCRIPT_STR
      "\n"
      "set(PROJECT_NAME ${PROJECT_NAME})\n"
      "\n"
      "set(${PROJECT_NAME}_TRIBITS_DIR ${${PROJECT_NAME}_TRIBITS_DIR})\n"
      "\n"
      "set(TEST_NAME ${TEST_NAME})\n"
      "\n"
      "set(NUM_CMNDS ${NUM_CMNDS})\n"
      "\n"
      "set(OVERALL_WORKING_DIRECTORY \"${PARSE_OVERALL_WORKING_DIRECTORY}\")\n"
      "\n"
      "set(SKIP_CLEAN_OVERALL_WORKING_DIRECTORY \"${PARSE_SKIP_CLEAN_OVERALL_WORKING_DIRECTORY}\")\n"
      "\n"
      "set(FAIL_FAST ${PARSE_FAIL_FAST})\n"
      "\n"
      "set(SHOW_START_END_DATE_TIME ${${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME})\n"
      "\n"
      "set(SHOW_MACHINE_LOAD ${${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST})\n"
      "\n"
      "set(CATEGORIES ${PARSE_CATEGORIES})\n"
      "\n"
      "set(PROCESSORS ${MAX_NUM_PROCESSORS_USED})\n"
      "\n"
      "set(TIMEOUT ${TIMEOUT_USED})\n"
      "\n"
      "#\n"
      "# Test invocation\n"
      "#\n"
      "\n"
      "include(\"${coreUtilsDir}/DriveAdvancedTest.cmake\")\n"
      "\n"
      "drive_advanced_test()\n"
      )
  
    if (TRIBITS_ADD_ADVANCED_TEST_UNITTEST)
      global_set(TRIBITS_ADD_ADVANCED_TEST_NUM_CMNDS ${NUM_CMNDS})
      # NOTE: This var only gets set if the advanced test gets added after
      # applying all of the various logic.  Therefore, unit tests should only
      # check this variable for being empty to determine that the test was not
      # added.
    endif()
  
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_var(TEST_SCRIPT_STR)
    endif()
  
    # Write the script file
  
    if (NOT TRIBITS_ADD_ADVANCED_TEST_SKIP_SCRIPT)
  
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("\nWriting file \"${TEST_SCRIPT_FILE}\" ...")
      endif()
  
      file( WRITE "${TEST_SCRIPT_FILE}"
        "${TEST_SCRIPT_STR}" )
  
    endif()


  else()

    global_set(TRIBITS_ADD_ADVANCED_TEST_NUM_CMNDS "")

  endif()

endfunction()
