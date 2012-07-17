# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

INCLUDE(TribitsAddTestHelpers)


#
# Add a test or a set of tests for a single executable.
#
#   TRIBITS_ADD_TEST(
#     <execName>
#     [NOEXEPREFIX]
#     [NOEXESUFFIX]
#     [NAME <testName> | NAME_POSTFIX <testNamePostfix>]
#     [DIRECTORY <directory>]
#     [ADD_DIR_TO_NAME]
#     [CREATE_WORKING_DIR]
#     [ARGS "<arg1> <arg2> ..." "<arg3> <arg4> ..." ...
#       | POSTFIX_AND_ARGS_0 <postfix> <arg1> <arg2> ...
#         POSTFIX_AND_ARGS_1 ... ]
#     [KEYWORDS <keyword1> <keyword2> ...]
#     [COMM [serial] [mpi]]
#     [NUM_MPI_PROCS <numProcs>]
#     [CATEGORIES <category1>  <category2> ...]
#     [HOST <host1> <host2> ...]
#     [XHOST <host1> <host2> ...]
#     [HOSTTYPE <hosttype1> <hosttype2> ...]
#     [XHOSTTYPE <hosttype1> <hosttype2> ...]
#     [STANDARD_PASS_OUTPUT
#       | PASS_REGULAR_EXPRESSION "<regex1>;<regex2>;..."]
#     [FAIL_REGULAR_EXPRESSION "<regex1>;<regex2>;..."]
#     [WILL_FAIL]
#     )
#  
# The arguments to the function are as followes:
#
#   <execName>
#
#     The mandatory root name of the executable that will be run to define the
#     test.  The full name of the executable is assumed to be
#     ${PACKAGE_NAME}_<execName>.exe and this executable is assumed to exist
#     in the current binary directory from where this function is called from
#     a CMakeLists.txt file.  This names is the default naming of executables
#     created by the function TRIBITS_ADD_EXECUTABLE(<execName> ...).
#     However, any arbitrary execuable program or script can be called by
#     setting NOEXEPREFIX and NOEXESUFFIX and you can also give <execName> as
#     an absolute path.
#
#   NOEXEPREFIX
#
#     If specified, then the prefix ${PACKAGE_NAME}_ is not assumed to be
#     appended to <execName>.
#
#   NOEXESUFFIX
#
#     If specified, then the postfix '.exe' is not assumed to be post-pended
#     to <execName>.
#
#   NAME <testName>
#
#     If specified, gives the root name of the test.  If not specified, then
#     <testName> is taked to be <execName>.  The actual test name will always
#     prefixed as ${PACKAGE_NAME}_<testName> and as added in the call to the
#     built-in CMake command ADD_TEST(...).  The main purpose of this argument
#     is to allow multiple tests to be defined for the same execurable.
#     Otherwise, you can not do this because CTest requires all test names to
#     be globally unique in a single project.
#
#   NAME_POSTFIX <testNamePostfix>
#
#     If specified, gives a postfix that will be added to the name of the
#     executable in order build the name of the test (appended as
#     '_<NAME_POSTFIX>').  If NAME argument is given it will be selected over
#     the NAME_POSTFIX.
#
#   DIRECTORY <directory>
#
#     If specified, then the executable is assumed to be in the directory
#     given by relative <directory>.
# 
#   ADD_DIR_TO_NAME
#
#     If specified the directory name that the test resides in will be added into
#     the name of the test after any package name is added and before the given
#     name of the test. the directory will have the package's base directory
#     stripped off so only the unique part of the test directory will be used.
#     All directory seperators will be changed into underscores.
#
#   CREATE_WORKING_DIR
#
#     If specified then a single temporary working directory will be created
#     based on the full name of the test and the test will be run in that
#     directory.  This allows tests to run in parallel even if they might
#     otherwise write the same files.  NOTE: If you use this option, then you
#     will need to adjust the path to your executable or other command.  Also,
#     input file locations will need to point down a directory.
#
#   ARGS "<arg1> <arg2> ..." "<arg3> <arg4> ..." ...
#
#     If specified, then a set of arguments can be passed in quotes.  If
#     multiple groups of arguments are passed in different quoted clusters of
#     arguments than a different test will be added for each set of arguments.
#     In this way, many different tests can be added for a single executable
#     in a single call to this function.  Each of these separate tests will be
#     named ${PACKAGE_NAME}_<testName>_xy where xy = 00, 01, 02, and so on.
#
#   POSTFIX_AND_ARGS_<IDX> <postfix> <arg1> <arg2> ...
#
#     If specified, gives a sequence of sets of test postfix names and
#     arguments lists for defining different tests.  For example, a set of
#     three different tests with argument lists can be specified as:
#     
#       POSTIFX_AND_ARGS_0 postfix1 --arg1 --arg2="dummy"
#       POSTIFX_AND_ARGS_1 postfix2  --arg2="fly"
#       POSTIFX_AND_ARGS_3 postfix3  --arg2="bags"
#
#     This will create three different test cases with the postfix names
#     'postfix1', 'postfix2', and 'postfix3'.  The indexes must be consecutive
#     starting a 0 and going up to (currently) 19.  The main advantages of
#     using these arguments instead of just 'ARGS' are that you can give
#     meaningful name to each test case and you can specify multiple arguments
#     without having to quote them and you can allow long argument lists to
#     span multiple lines.
#
#   KEYWORDS <keyword1> <keyword2> ...
#
#     If specified, gives a list of keywords added to a test.  These keywords
#     can then be used to select tests to be run with 'ctest'.
#
#   COMM [serial] [mpi]
#
#     If specified, selects if the test will be added in serial and/or MPI
#     mode.  If the COMM argument is missing, the test will be added in both
#     serial and MPI builds of the code.
#
#   NUM_MPI_PROCS <numProcs>
#
#     If specified, gives the number of processes that the test will be
#     defined to run on and can also result in the test being exluded all
#     together based on comparison to MPI_EXEC_MAX_NUMPROCS.  *** ToDo: Finish
#     specification of this arugment! ***
#
#   HOST <host1> <host2> ...
#
#     If specified, gives a list of hostnames where the test will be included.
#     The current hostname is determined by the built-in CMake command
#     SITE_NAME(...).  On Linux/Unix systems, this is typically the value
#     returned by 'uname -n'.
#
#   XHOST <host1> <host2> ...
#
#     If specified, gives a list of hostnames where the test will *not* be
#     included.  This check is performed after the check for the hostnames in
#     the 'HOST' list if it should exist.  Therefore, this list exclusion list
#     overrides the 'HOST' inclusion list.
#
#   CATEGORIES <category1> <category2> ...
#
#     If specified, gives the specific categories of the test.  Valid test
#     categories include BASIC, CONTINUOUS, NIGHTLY, and PERFORMANCE.  Other
#     test categories will be added as needed.  By default, the category is
#     BASIC.  When the CATEGORIES is BASIC it will match
#     ${PROJECT_NAME}_TEST_CATEGORIES equal to CONTINUOUS and NIGHTLY.  When
#     the CATEGORIES is CONTINUOUS it will match
#     ${PROJECT_NAME}_TEST_CATEGORIES equal to CONTINUOUS and NIGHTLY.  When
#     the CATEGORIES is PERFORMANCE it will match
#     ${PROJECT_NAME}_TEST_CATEGORIES=PERFORMANCE ony.
#
#   HOSTTYPE <hosttype1> <hosttype2> ...
#
#     If specified, gives the names of the host system type (given by
#     ${CMAKE_HOST_SYSTEM_NAME}) to include the test.  Typical host system
#     type names include 'Linux', 'Darwain' etc.
#
#   XHOSTTYPE <name1> <name2> ...
#
#     If specified, gives the names of the host system type to *not* include
#     the test.  This check is performed after the check for the host system
#     names in the 'HOSTTYPE' list if it should exist.  Therefore, this list
#     exclusion list overrides the 'HOSTTYPE' inclusion list.
#
#   STANDARD_PASS_OUTPUT
#
#     If specified, then the standard test output "End Result: TEST PASSED" is
#     greped for to determine success.  This is needed for MPI tests on some
#     platforms since the return value is unreliable.
#
#   PASS_REGULAR_EXPRESSION "<regex1>;<regex2>;..." 
#
#     If specified, then a test will be assumed to pass only if one of the
#     regular expressions <regex1>, <regex2> etc. match the output.
#     Otherwise, the test will fail.
#
#   FAIL_REGULAR_EXPRESSION "<regex1>;<regex2>;..."
#
#     If specified, then a test will be assumed to fail if one of the regular
#     expressions <regex1>, <regex2> etc. match the output.  Otherwise, the
#     test will pass.
#
#   WILL_FAIL
#
#     If passed in, then the pass/fail criteria will be inverted.
#
# NOTES:
#
# 1) The test can be disabled by setting the variable
# ${PACKAGE_NAME}_${TEST_NAME}_DISABLE=TRUE (perhaps in the cache).  This allows
# tests to be disable on a case-by-case basis.  Here, TEST_NAME is ${EXE_NAME}
# plus any other postfixes.  This is the name that shows up in 'ctest -N' when
# running the test.
#

FUNCTION(TRIBITS_ADD_TEST EXE_NAME)

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRIBITS_ADD_TEST: ${EXE_NAME} ${ARGN}")
  ENDIF()

  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT "")
   
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
     "DIRECTORY;KEYWORDS;COMM;NUM_MPI_PROCS;ARGS;${POSTFIX_AND_ARGS_LIST};NAME;NAME_POSTFIX;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;PASS_REGULAR_EXPRESSION;FAIL_REGULAR_EXPRESSION;TIMEOUT"
     #options
     "NOEXEPREFIX;NOEXESUFFIX;STANDARD_PASS_OUTPUT;WILL_FAIL;ADD_DIR_TO_NAME;CREATE_WORKING_DIR"
     ${ARGN}
     )
  # NOTE: The TIMEOUT argument is not documented on purpose.  I don't want to
  # advertise it!

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

  TRIBITS_ADD_TEST_GET_EXE_BINARY_NAME( "${EXE_NAME}"
    ${PARSE_NOEXEPREFIX} ${PARSE_NOEXESUFFIX} ${PARSE_ADD_DIR_TO_NAME} EXE_BINARY_NAME )
  
  #If requested create a modifier for the name that will be inserted between the package name 
  #and the given name or exe_name for the test
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
        "${EXECUTABLE_PATH}"  "${NUM_PROCS_USED}" ${PARSE_CREATE_WORKING_DIR}
        ${INARGS} )
  
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
        ${INARGS} )

    ENDFOREACH()

  ENDIF()
  
ENDFUNCTION()
