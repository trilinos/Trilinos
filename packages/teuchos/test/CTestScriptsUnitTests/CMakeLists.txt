
#
# The tests defined here are used to test the dynamic behavior of the
# CMake/CTest functions PACKAGE_ADD_TEST(...) and
# PACKAGE_ADD_ADVANCED_TEST(...).  The reason that these tests are defined
# here in Teuchos is so that MPI handling can be taken care of automatically.
# These tests run extremely fast so they should not be a drain on testing of
# Teuchos or the rest of Trilinos.
#


INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddAdvancedTest)

PACKAGE_ADD_EXECUTABLE(
  GenericDriver
  SOURCES GenericDriver.cpp
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_0
  TEST_0 EXEC GenericDriver
  )

PACKAGE_ADD_TEST( GenericDriver
  NAME CTestScripts_simple_pass_regex
  ARGS "--print-msg=\"This is a crazy test!\" --return-val=1"
  PASS_REGULAR_EXPRESSION "[^=]This is a crazy test"
  )

PACKAGE_ADD_TEST( GenericDriver
  NAME CTestScripts_simple_posfix_and_args_1
  POSTFIX_AND_ARGS_0 test1 --print-msg="This is a crazy test!" --return-val=1
  PASS_REGULAR_EXPRESSION "[^=]This is a crazy test"
  )

PACKAGE_ADD_TEST( GenericDriver
  NAME CTestScripts_simple_posfix_and_args_2
  POSTFIX_AND_ARGS_0 test1 --print-msg="This is a crazy test!"
  POSTFIX_AND_ARGS_1 test2 --print-msg="This is a silly test!"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_print_msg
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test!" --return-val=1
  FINAL_PASS_REGULAR_EXPRESSION "[^=]This is a crazy test"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_print_msg_no_echo_output
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test!"
    NO_ECHO_OUTPUT
  FINAL_FAIL_REGULAR_EXPRESSION "[^=]This is a crazy test"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_return_code
  TEST_0 EXEC GenericDriver ARGS --return-val=1
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: Return code = 1"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_print_msg_out_file
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test!" --return-val=5
    OUTPUT_FILE cmnd_1_args_1_print_msg_out_file_outputFile.out
    PASS_REGULAR_EXPRESSION "This is a crazy test"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_any_test_0
  TEST_0 EXEC GenericDriver ARGS --return-val=7
    PASS_ANY
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: Pass criteria = Pass Any"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_any_overall
  TEST_0 EXEC GenericDriver ARGS --return-val=7
    PASS_ANY
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_standard_pass_output
  TEST_0 EXEC GenericDriver ARGS --return-val=7
    STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_regular_expression_all
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test\n\nThis is not the best test" --return-val=5
    OUTPUT_FILE cmnd_1_args_1_pass_regular_expression_all_outputFile.out
    PASS_REGULAR_EXPRESSION_ALL "This is a crazy test" "This is not the best test"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_regular_expression_all_compare_1
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test\n\nThis is not the best test" --return-val=5
    OUTPUT_FILE cmnd_1_args_1_pass_regular_expression_all_compare_1_outputFile.out
    PASS_REGULAR_EXPRESSION_ALL "This is a crazy test" "This is not the best test"
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: Pass criteria = Match REGEX {This is a crazy test} [^a]PASSED[^a]"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_pass_regular_expression_all_compare_2
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test\n\nThis is not the best test" --return-val=5
    OUTPUT_FILE cmnd_1_args_1_pass_regular_expression_all_compare_2_outputFile.out
    PASS_REGULAR_EXPRESSION_ALL "This is a crazy test" "This is not the best test"
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: Pass criteria = Match REGEX {This is not the best test} [^a]PASSED[^a]"
  )

# 2009/04/06: rabartl: Above, I can't use '[' or ']' in the regular expression
# because it is interpreted as a char list.  Therefore, I leave it as not 'a'
# and that is good enough for me.

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_print_msg_out_file_no_echo_output
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test!" --return-val=5
    OUTPUT_FILE cmnd_1_args_1_print_msg_out_file_no_echo_output_outputFile.out
    NO_ECHO_OUTPUT
  FINAL_FAIL_REGULAR_EXPRESSION "[^=]This is a crazy test"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_message
  TEST_0 EXEC GenericDriver ARGS --return-val=1 MESSAGE "This is a crazy test!"
  FINAL_PASS_REGULAR_EXPRESSION "This is a crazy test"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_1_args_1_working_directory
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test!"
    WORKING_DIRECTORY TEST_NAME
    OUTPUT_FILE cmnd_1_args_1_working_directory_console.out
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_args_1_basic
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is the first test" --return-val=3
    PASS_REGULAR_EXPRESSION "first test"
  TEST_1 EXEC GenericDriver ARGS --print-msg="This is the second test" --return-val=4
    PASS_REGULAR_EXPRESSION "second test"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_args_1_check_second_print
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is the first test" --return-val=3
    PASS_REGULAR_EXPRESSION "first test"
  TEST_1 EXEC GenericDriver ARGS --print-msg="This is the second test" --return-val=4
    PASS_REGULAR_EXPRESSION "second test"
  FINAL_PASS_REGULAR_EXPRESSION "second test"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_fail_fast
  TEST_0 EXEC GenericDriver ARGS --print-msg="This first test will fail" --return-val=3
  TEST_1 EXEC GenericDriver ARGS --print-msg="This is the second test will never be run"
  FAIL_FAST
  FINAL_PASS_REGULAR_EXPRESSION "TEST_0: FAIL FAST, SKIPPING REST OF TEST CASES!"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_no_fail_fast_print_msg
  TEST_0 EXEC GenericDriver ARGS --print-msg="This first test will fail" --return-val=3
  TEST_1 EXEC GenericDriver ARGS --print-msg="This is the second test should be run"
  FINAL_PASS_REGULAR_EXPRESSION "This is the second test should be run"
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_2_no_fail_fast_final_fail
  TEST_0 EXEC GenericDriver ARGS --print-msg="This first test will fail" --return-val=3
  TEST_1 EXEC GenericDriver ARGS --print-msg="This is the second test should be run"
  FINAL_PASS_REGULAR_EXPRESSION "OVERALL FINAL RESULT: TEST FAILED"
  )

IF (WIN23 AND NOT CYGWIN)
  SET(DIFF_COMMAND_STUFF "${CMAKE_COMMAND}" ARGS -E compare_files)
ELSE()
  SET(DIFF_COMMAND_STUFF diff ARGS)
ENDIF()

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_3_args_1_overall_working_directory
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test!"
    OUTPUT_FILE cmnd_3_args_1_overall_working_directory_console.out
  TEST_1 CMND "${CMAKE_COMMAND}" ARGS -E copy cmnd_3_args_1_overall_working_directory_console.out cmnd_3_args_1_overall_working_directory_console.gold.out
  TEST_2 CMND ${DIFF_COMMAND_STUFF} cmnd_3_args_1_overall_working_directory_console.out cmnd_3_args_1_overall_working_directory_console.gold.out
  )

PACKAGE_ADD_ADVANCED_TEST(
  CTestScripts_cmnd_3_args_1_overall_working_directory_working_directory
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 EXEC GenericDriver ARGS --print-msg="This is a crazy test!"
    OUTPUT_FILE cmnd_3_args_1_overall_working_directory_working_directory_console.out
    WORKING_DIRECTORY GenericDriver
  TEST_1 CMND "${CMAKE_COMMAND}" ARGS -E copy ../GenericDriver/cmnd_3_args_1_overall_working_directory_working_directory_console.out cmnd_3_args_1_overall_working_directory_working_directory_console.gold.out
    WORKING_DIRECTORY Gold
  TEST_2 CMND ${DIFF_COMMAND_STUFF} GenericDriver/cmnd_3_args_1_overall_working_directory_working_directory_console.out Gold/cmnd_3_args_1_overall_working_directory_working_directory_console.gold.out
  )

# Demonstrate processes never dying

IF (PERL_EXECUTABLE AND Teuchos_ENABLE_RUNS_FOREVER_TEST)
  
  PACKAGE_ADD_EXECUTABLE(
    RunsForever
    SOURCES RunsForever.cpp
    )
  
  PACKAGE_ADD_ADVANCED_TEST(
    CTestScripts_RunsForever
    TEST_0 CMND ${PERL_EXECUTABLE} ARGS ${CMAKE_CURRENT_SOURCE_DIR}/runs_forever.pl
    )

ENDIF()
