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


#
# The tests defined here are used to test the dynamic behavior of the
# CMake/CTest functions TRIBITS_ADD_TEST(...) and
# TRIBITS_ADD_ADVANCED_TEST(...).  The reason that these tests are defined
# here in Teuchos is so that MPI handling can be taken care of automatically.
# These tests run extremely fast so they should not be a drain on testing of
# Teuchos or the rest of Trilinos.
#


INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddAdvancedTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_MPI)
IF (${PACKAGE_NAME}_ENABLE_MPI)
  ADD_DEFINITIONS(-DHAVE_MPI)
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  GenericDriver
  SOURCES GenericDriver.cpp
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_0
  TEST_0 EXEC GenericDriver
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST( GenericDriver
  NAME CTestScripts_simple_pass_regex
  ARGS "\"This is a crazy test!\" 1"
  PASS_REGULAR_EXPRESSION "[^=]This is a crazy test"
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST( GenericDriver
  NAME CTestScripts_simple_posfix_and_args_1
  POSTFIX_AND_ARGS_0 test1 "This is a crazy test!" 1
  PASS_REGULAR_EXPRESSION "[^=]This is a crazy test"
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST( GenericDriver
  NAME CTestScripts_simple_posfix_and_args_2
  POSTFIX_AND_ARGS_0 test1 "This is a crazy test!" 0
  POSTFIX_AND_ARGS_1 test2 "This is a silly test!" 0
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST( GenericDriver
  NAME CTestScripts_will_fail
  POSTFIX_AND_ARGS_0 test1 "This test is supposed to fail!" 1
  WILL_FAIL
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_print_msg
  TEST_0 EXEC GenericDriver ARGS "This is a crazy test!" 1
  FINAL_PASS_REGULAR_EXPRESSION "\"This is a crazy test"
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_print_msg_no_echo_output
  TEST_0 EXEC GenericDriver ARGS "This is a crazy test!" 0
    NO_ECHO_OUTPUT
  FINAL_FAIL_REGULAR_EXPRESSION "\"This is a crazy test\""
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_return_code
  TEST_0 EXEC GenericDriver ARGS "Dummy msg" 1
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: Return code = 1"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_print_msg_out_file
  TEST_0 EXEC GenericDriver ARGS "This is a crazy test!" 5
    OUTPUT_FILE cmnd_1_args_1_print_msg_out_file_outputFile.out
    PASS_REGULAR_EXPRESSION "This is a crazy test"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_any_test_0
  TEST_0 EXEC GenericDriver ARGS "Dummy msg" 7
    PASS_ANY
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: Pass criteria = Pass Any"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_any_overall
  TEST_0 EXEC GenericDriver ARGS "Dummy msg" 7
    PASS_ANY
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_standard_pass_output
  TEST_0 EXEC GenericDriver ARGS "Dummy msg" 7
    STANDARD_PASS_OUTPUT
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_2_stderr_pass_regular_expression_all
  TEST_0 EXEC GenericDriver
    ARGS "This is a crazy test\n\nThis is not the best test" 5 "This goes to stderr"
    OUTPUT_FILE cmnd_1_args_1_pass_regular_expression_all_outputFile.out
    PASS_REGULAR_EXPRESSION_ALL "This is a crazy test" "This is not the best test" "This goes to stderr"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_regular_expression_all_compare_1
  TEST_0 EXEC GenericDriver
    ARGS "This is a crazy test\n\nThis is not the best test" 5
    OUTPUT_FILE cmnd_1_args_1_pass_regular_expression_all_compare_1_outputFile.out
    PASS_REGULAR_EXPRESSION_ALL "This is a crazy test" "This is not the best test"
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: Pass criteria = Match REGEX {This is a crazy test} [^a]PASSED[^a]"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_regular_expression_all_compare_2
  TEST_0 EXEC GenericDriver
    ARGS "This is a crazy test\n\nThis is not the best test" 5
    OUTPUT_FILE cmnd_1_args_1_pass_regular_expression_all_compare_2_outputFile.out
    PASS_REGULAR_EXPRESSION_ALL "This is a crazy test" "This is not the best test"
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: Pass criteria = Match REGEX {This is not the best test} [^a]PASSED[^a]"
  OVERALL_NUM_MPI_PROCS 1
  )

# 2009/04/06: rabartl: Above, I can't use '[' or ']' in the regular expression
# because it is interpreted as a char list.  Therefore, I leave it as not 'a'
# and that is good enough for me.

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_print_msg_out_file_no_echo_output
  TEST_0 EXEC GenericDriver
    ARGS "This is a crazy test!" 5
    OUTPUT_FILE cmnd_1_args_1_print_msg_out_file_no_echo_output_outputFile.out
    NO_ECHO_OUTPUT
  FINAL_FAIL_REGULAR_EXPRESSION "^This is a crazy test"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_message
  TEST_0 EXEC GenericDriver
    ARGS "Dummy msg" 1
    MESSAGE "This is a crazy test!"
  FINAL_PASS_REGULAR_EXPRESSION "This is a crazy test"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_working_directory
  TEST_0 EXEC GenericDriver
    ARGS "This is a crazy test!" 0
    WORKING_DIRECTORY TEST_NAME
    OUTPUT_FILE cmnd_1_args_1_working_directory_console.out
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_args_1_basic
  TEST_0 EXEC GenericDriver
    ARGS "This is the first test" 3
    PASS_REGULAR_EXPRESSION "first test"
  TEST_1 EXEC GenericDriver
    ARGS "This is the second test" 4
    PASS_REGULAR_EXPRESSION "second test"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_args_1_check_second_print
  TEST_0 EXEC GenericDriver
    ARGS "This is the first test" 3
    PASS_REGULAR_EXPRESSION "first test"
  TEST_1 EXEC GenericDriver
    ARGS "This is the second test" 4
    PASS_REGULAR_EXPRESSION "second test"
  FINAL_PASS_REGULAR_EXPRESSION "second test"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_fail_fast
  TEST_0 EXEC GenericDriver
    ARGS "This first test will fail" 3
  TEST_1 EXEC GenericDriver
    ARGS "This is the second test will never be run" 0
  FAIL_FAST
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: FAIL FAST, SKIPPING REST OF TEST CASES!"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_no_fail_fast_print_msg
  TEST_0 EXEC GenericDriver
    ARGS "This first test will fail" 3
  TEST_1 EXEC GenericDriver
    ARGS "This is the second test should be run" 0
  FINAL_PASS_REGULAR_EXPRESSION "This is the second test should be run"
  OVERALL_NUM_MPI_PROCS 1
  )

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_no_fail_fast_final_fail
  TEST_0 EXEC GenericDriver
    ARGS "This first test will fail" 3
  TEST_1 EXEC GenericDriver
    ARGS "This is the second test should be run" 0
  FINAL_PASS_REGULAR_EXPRESSION "OVERALL FINAL RESULT: TEST FAILED"
  OVERALL_NUM_MPI_PROCS 1
  )

IF (CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
  # cmake -E compare_files only gives return != 0 if the files don't compare.
  # It gives no diff output.
  SET(DIFF_COMMAND_STUFF "${CMAKE_COMMAND}" ARGS -E compare_files)
ELSE()
  # Diff gives output of the diff so it is better to use if we have it.
  # Otherwise, I need to write a portable diff (using a CMake -P script).
  SET(DIFF_COMMAND_STUFF diff ARGS)
ENDIF()

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_3_args_1_overall_working_directory
  MESSAGE "NOTE: This test does not check the order of STDOUT and STDERR, only that both exist in the output file.  It seems that EXECUTE_PROCESS() and MPI do not behave well together w.r.t. to STDOUT and STDERR outputting.  The fault may be with the MPI implementation.  Also, even in serial it would seem that the order of STDOUT and STDERR gets messed up when the machine is loaded."
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 EXEC GenericDriver
    ARGS "This is a crazy test!" 0 "This is a bad crazy test!"
    OUTPUT_FILE cmnd_3_args_1_overall_working_directory_console.out
  TEST_1 CMND cat
    ARGS cmnd_3_args_1_overall_working_directory_console.out
    PASS_REGULAR_EXPRESSION_ALL
       "This is a crazy test"
       "This is a bad crazy test"
  OVERALL_NUM_MPI_PROCS 1
  XHOSTTYPE Windows
  )
# NOTE: The above test does not run on Windows because it relies on the 'cat'
# command.  The form of this test is needed in order to check that both STDOUT
# and STDERR output are both in the output file but we can't rely on the
# order.  It seems that EXECUTE_PROCESS() and mpiexec don't work together to
# produce ordered STDOUT and STDERR very well.  Also, even in serial, it seems
# that when the machine is loaded that the order of STDOUT and STDERR can get
# messed up.  It is not clear if this is a problem with EXECUTE_PROCESS() or
# not.

TRIBITS_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_3_args_1_overall_working_directory_working_directory
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 EXEC GenericDriver
    ARGS "This is another crazy test!" 0
    OUTPUT_FILE cmnd_3_args_1_overall_working_directory_working_directory_console.out
    WORKING_DIRECTORY GenericDriver
  TEST_1 CMND "${CMAKE_COMMAND}"
    ARGS -E copy
      ../GenericDriver/cmnd_3_args_1_overall_working_directory_working_directory_console.out
      cmnd_3_args_1_overall_working_directory_working_directory_console.gold.out
    WORKING_DIRECTORY Gold
  TEST_2 CMND ${DIFF_COMMAND_STUFF}  # Contains the ARGS part
    GenericDriver/cmnd_3_args_1_overall_working_directory_working_directory_console.out
    Gold/cmnd_3_args_1_overall_working_directory_working_directory_console.gold.out
  OVERALL_NUM_MPI_PROCS 1
  )

# Demonstrate processes never dying

IF (PERL_EXECUTABLE AND TrilinosFramework_ENABLE_RUNS_FOREVER_TEST)
  
  TRIBITS_ADD_EXECUTABLE(
    RunsForever
    SOURCES RunsForever.cpp
    )
  
  TRIBITS_ADD_ADVANCED_TEST(
    CTestScripts_RunsForever
    TEST_0 CMND ${PERL_EXECUTABLE} ARGS ${CMAKE_CURRENT_SOURCE_DIR}/runs_forever.pl
    OVERALL_NUM_MPI_PROCS 1
    )

ENDIF()

# ENVIRONMENT

TRIBITS_ADD_TEST( "${CMAKE_COMMAND}" NOEXEPREFIX NOEXESUFFIX
  NAME CTestScripts_TAT_pass_env
  ARGS "-P ${CMAKE_CURRENT_SOURCE_DIR}/EchoEnvVarForTest.cmake"
  NUM_MPI_PROCS 1
  ENVIRONMENT PASS_IN_FROM_ENV_VAR=ThisWasPassedCorrectly
  PASS_REGULAR_EXPRESSION "PASS_IN_FROM_ENV_VAR: ThisWasPassedCorrectly"
  )

TRIBITS_ADD_ADVANCED_TEST( CTestScripts_TAAT_pass_env
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS -P ${CMAKE_CURRENT_SOURCE_DIR}/EchoEnvVarForTest.cmake
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "PASS_IN_FROM_ENV_VAR: ThisWasPassedCorrectly"
  ENVIRONMENT PASS_IN_FROM_ENV_VAR=ThisWasPassedCorrectly
  )
