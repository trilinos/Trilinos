
INCLUDE(PackageAddAdvancedTestHelpers)

INCLUDE(AppendStringVar)
INCLUDE(PrintVar)


#
# Function that creates an advanced test defined using one or more executable
# commands that is run as a separate CMake script.
#
# This function allows you to add a single CTest test as a single unit that is
# actually a sequence of one or more separate commands strung together in some
# way to define the final pass/fail.
#
# An advanced test is defined as:
#
#   PACKAGE_ADD_ADVANCED_TEST(
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
#     [CATEGORIES <category1> <category2> ...]
#     [HOST <host1> <host2> ...]
#     [XHOST <host1> <host2> ...]
#     [FINAL_PASS_REGULAR_EXPRESSION <regex> | FINAL_FAIL_REGULAR_EXPRESSION <regex>]
#     )
#
# Each and every atomic test or command needs to pass (as defined below) in
# order for the overall test to pass
#
# Each atomic test case is either a package-built executable or just a basic
# command.  An atomic test command takes the form:
#
#   TEST_<i>
#      EXEC <execTarget> [NOEXEPREFIX] [NOEXESUFFIX] [ADD_DIR_TO_NAME]
#         | CMND <cmndExec>
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
# ToDO: Add documnetation for [X]HOST[TYPE]
#
# Some overall arguments are:
#
#   <testName>
#
#     The name of the test (which will have ${PACKAGE_NAME}_
#     appended) that will be used to name the output CMake script file as well as
#     the CTest test name passed into ADD_TEST(...).
#
#   TEST_<i> (EXEC <execTarget0> | CMND <cmndExec0>) ...
#
#     Defines test command <i>.  Each of these test commands must be in
#     sequential order.  The details for each atomic test are given below.
#
#   OVERALL_WORKING_DIRECTORY <overallWorkingDir>
#
#     If specified, then the working directory <overallWorkingDir> will be
#     created and all of the test commands by default will be run from within
#     this directory.  If the value <overallWorkingDir> = TEST_NAME is given,
#     then the working directory will be given the name
#     ${PACKAGE_NAME}_<testName>.  If the directory <overallWorkingDir> exists
#     before the test runs, it will be deleted and created again.  Therefore,
#     if you want to preserve the contents of this directory between test runs
#     you need to copy it somewhere else.
#
#   KEYWORDS <keyword1> <keyword2> ...
#
#     If specified, gives a list of keywords added to a test.  These keywords
#     can then be used to select tests to be run with 'ctest'.
#
#   FAIL_FAST
#
#     If specified, then the remaining test commands will be aborted when any
#     test command fails.  Otherwise, all of the test cases will be run.
#
#   COMM [serial] [mpi]
#
#     If specified, selects if the test will be added in serial and/or MPI
#     mode.  If the COMM argument is missing, the test will be added in both
#     serial and MPI builds of the code.  See the COMM argument in the script
#     PACKAGE_ADD_TEST(...) for more details.
#
#   OVERALL_NUM_MPI_PROCS <overallNumProcs>
#
#     If specified, gives the default number of processes that each executable
#     command run on and can also result in the test being exluded all
#     together based on comparison to MPI_EXEC_MAX_NUMPROCS.  See the COMM
#     argument in the script PACKAGE_ADD_TEST(...) for more details.
#
# Each test command is either package-built test executable or some general
# command executable and is defined as either:
#
#   EXEC <execTarget>
#
#     If specified, then <execTarget> gives the the name of an executable
#     target that will be run as the command.  The value <execTarget> is same
#     string that was passed in as the first argument to
#     PACKAGE_ADD_EXECUTABLE( <execTarget>...) used to define the executable.
#     If this is an MPI build, then <execTarget> will be run with MPI using
#     NUM_MPI_PROCS <numProcs> or OVERALL_NUM_MPI_PROCS <overallNumProcs> (if
#     NUM_MPI_PROCS is not set for this test case).  If the number of maximum
#     MPI processes allowed is less than this number of MPI processes, then
#     the test will *not* be run.  If NOEXEPREFIX is specified, then
#     ${PACKAGE_NAME}_ will not be added the the beginning.  If NOEXESUFFIX is
#     specified, then '.exe' will not be added to the end.  Note that EXEC
#     <execTarget> is basically equivalent to CMND <cmndExec> when NOEXEPREFIX
#     and NOEXESUFFIX are specified.  In this case, you can pass in
#     <execTarget> to any command you would like and it will get run with MPI
#     in MPI mode just link any other command.
#
#   CMND <cmndExec>
#
#     If specified, then <cmndExec> gives the executable for a command to be
#     run.  In this case, MPI will never be used to run the executable even
#     when configured in MPI mode (i.e. TPL_ENABLE_MPI=ON).
#
# By defualt, the output (stdout/stderr) for each test command is captured and
# is then echoed to stdout for the overall test.  This is done in order to be
# able to grep the result to determine pass/fail.
#
# Other miscellaneous arguments include:
#
#   MESSAGE "<message>"
#
#     If specified, then the string in <message> will be print before this
#     test command is run.
#
#   WORKING_DIRECTORY <workingDir>
#
#     If specified, then the working directory <workingDir> will be created
#     and the test will be run from within this directory.  If the value
#     <workingDir> = TEST_NAME is given, then the working directory will be
#     given the name ${PACKAGE_NAME}_<testName>.  If the directory
#     <workingDir> exists before the test runs, it will be deleted and created
#     again.  Therefore, if you want to preserve the contents of this
#     directory between test runs you need to copy it somewhere else.
#
#   NUM_MPI_PROCS <numProcs>
#
#     If specified, then <numProcs> is the number of processors used for MPI
#     executables.  If not specified, this will default to <overallNumProcs>
#     from OVERALL_NUM_MPI_PROCS <overallNumProcs>.
#
#   OUTPUT_FILE <outputFile>
#
#     If specified, then stdout and stderr will be sent to <outputFile>.
#
#   NO_ECHO_OUTPUT
#
#     If specified, then the output for the test command will not be echoed to
#     the output for the entire test command.
#
# By default, an atomic test line is assumed to pass if the executable returns
# a non-zero value.  However, a test case can also be defined to pass based
# on:
#
#   PASS_ANY
#
#     If specified, the test command 'i' will be assumed to pass reguardless
#     of the return value or any other output.  This would be used when a
#     command that is to follow will determine pass or fail based on output
#     from this command in some way.
#
#   PASS_REGULAR_EXPRESSION "<regex>"
#
#     If specified, the test command 'i' will be assumed to pass if it matches
#     the given regular expression.  Otherwise, it is assumed to fail.
#
#   PASS_REGULAR_EXPRESSION_ALL "<regex1>" "<regex2>" ... "<regexn>"
#
#     If specified, the test command 'i' will be assumed to pas if the output
#     matches all of the provided regular expressions.
#
#   FAIL_REGULAR_EXPRESSION "<regex>"
#
#     If specified, the test command 'i' will be assumed to fail if it matches
#     the given regular expression.  Otherwise, it is assumed to pass.
#
#   STANDARD_PASS_OUTPUT
#
#     If specified, the test command 'i' will be assumed to pass if the string
#     expression "Final Result: PASSED" is found in the ouptut for the test.
#
# By default, the overall test will be assumed to pass if it prints:
#
#  "OVERALL FINAL RESULT: TEST PASSED"
#
# However, this can be changed by setting one of the following optional arguments:
#
#   FINAL_PASS_REGULAR_EXPRESSION <regex>
#
#     If specified, the test will be assumed to pass if the output matches
#     <regex>.  Otherwise, it will be assumed to fail.
#
#   FINAL_FAIL_REGULAR_EXPRESSION <regex>
#
#     If specified, the test will be assumed to fail if the output matches
#     <regex>.  Otherwise, it will be assumed to fail.
#

FUNCTION(PACKAGE_ADD_ADVANCED_TEST TEST_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_ADD_ADVANCED_TEST: ${TEST_NAME_IN}\n")
  ENDIF()

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
  SET(MAX_NUM_TEST_CMND_IDX ${PACKAGE_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX})

  SET(TEST_IDX_LIST "")
  FOREACH( TEST_CMND_IDX RANGE ${MAX_NUM_TEST_CMND_IDX})
    LIST( APPEND TEST_IDX_LIST TEST_${TEST_CMND_IDX} )
  ENDFOREACH()
  #PRINT_VAR(TEST_IDX_LIST)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "${TEST_IDX_LIST};OVERALL_WORKING_DIRECTORY;KEYWORDS;COMM;OVERALL_NUM_MPI_PROCS;FINAL_PASS_REGULAR_EXPRESSION;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;FINAL_FAIL_REGULAR_EXPRESSION"
     #options
     "FAIL_FAST"
     ${ARGN}
     )
  
  #
  # B) Add or don't add tests based on a number of criteria
  #

  SET(ADD_THE_TEST FALSE)
  PACKAGE_ADD_TEST_PROCESS_CATEGORIES(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  SET(ADD_THE_TEST FALSE)
  PACKAGE_ADD_TEST_PROCESS_HOST_HOSTTYPE(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  #
  # C) Determine if we will add the serial or MPI tests based on input COMM
  # and TPL_ENABLE_MPI
  #

  PACKAGE_PROCESS_COMM_ARGS(ADD_SERIAL_TEST  ADD_MPI_TEST  ${PARSE_COMM})
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
       "EXEC;CMND;ARGS;MESSAGE;WORKING_DIRECTORY;OUTPUT_FILE;NUM_MPI_PROCS;PASS_REGULAR_EXPRESSION_ALL;FAIL_REGULAR_EXPRESSION;PASS_REGULAR_EXPRESSION"
       #options
       "NOEXEPREFIX;NOEXESUFFIX;NO_ECHO_OUTPUT;PASS_ANY;STANDARD_PASS_OUTPUT;ADD_DIR_TO_NAME"
       ${PARSE_TEST_${TEST_CMND_IDX}}
       )

    # Write the command

    SET(ARGS_STR ${PARSE_ARGS})
    #PRINT_VAR(ARGS_STR)
    #IF (PARSE_ARGS)
    #  JOIN_EXEC_PROCESS_SET_ARGS( ARGS_STR ${PARSE_ARGS} )
    #ENDIF()

    IF (PARSE_EXEC)

      LIST( LENGTH PARSE_EXEC PARSE_EXEC_LEN )
      IF (NOT PARSE_EXEC_LEN EQUAL 1)
        MESSAGE(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} EXEC = '${PARSE_EXEC}'"
          " must be a single name.  To add arguments use ARGS <arg1> <arg2> ...." )
      ENDIF()

      PACKAGE_ADD_TEST_GET_EXE_BINARY_NAME( "${PARSE_EXEC}"
        ${PARSE_NOEXEPREFIX} ${PARSE_NOEXESUFFIX} ${PARSE_ADD_DIR_TO_NAME} EXE_BINARY_NAME )

      IF (IS_ABSOLUTE ${EXE_BINARY_NAME})
        SET(EXECUTABLE_PATH "${EXE_BINARY_NAME}")
      ELSE()
        SET(EXECUTABLE_PATH "${CMAKE_CURRENT_BINARY_DIR}/${EXE_BINARY_NAME}")
      ENDIF()

      IF (NOT PARSE_NUM_MPI_PROCS)
        SET(PARSE_NUM_MPI_PROCS ${PARSE_OVERALL_NUM_MPI_PROCS})
      ENDIF()
      PACKAGE_ADD_TEST_GET_NUM_PROCS_USED("${PARSE_NUM_MPI_PROCS}" NUM_PROCS_USED)
      IF (NUM_PROCS_USED LESS 0)
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE(STATUS "Skipping adding test because NUM_MPI_PROCS ${PARSE_NUM_MPI_PROCS}"
            " is out of range.")
        ENDIF()
        SET(ADD_THE_TEST FALSE)
      ENDIF()

      PACKAGE_ADD_TEST_GET_TEST_CMND_ARRAY( TEST_CMND_ARRAY
        "${EXECUTABLE_PATH}" "${NUM_PROCS_USED}" ${ARGS_STR} )
      #PRINT_VAR(TEST_CMND_ARRAY)

    ELSEIF (PARSE_CMND)

      LIST( LENGTH PARSE_CMND PARSE_CMND_LEN )
      IF (NOT PARSE_CMND_LEN EQUAL 1)
        MESSAGE(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} CMND = '${PARSE_CMND}'"
          " must be a single command.  To add arguments use ARGS <arg1> <arg2> ...." )
      ENDIF()

      #This allows us to check if a test is requesting more than MPI_EXEC_MAX_NUMPROCS
      PACKAGE_ADD_TEST_GET_NUM_PROCS_USED("${PARSE_OVERALL_NUM_MPI_PROCS}" NUM_PROCS_USED)

      IF (NUM_PROCS_USED LESS 0)
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE(STATUS "Skipping adding test because OVERALL_NUM_MPI_PROCS ${PARSE_OVERALL_NUM_MPI_PROCS}"
            " is out of range.")
        ENDIF()
        SET(ADD_THE_TEST FALSE)
      ENDIF()

      SET( TEST_CMND_ARRAY ${PARSE_CMND} ${ARGS_STR} )
    
    ELSE()

      MESSAGE( FATAL_ERROR
        "Must have EXEC or CMND for TEST_${TEST_CMND_IDX}" )

    ENDIF()

    JOIN_EXEC_PROCESS_SET_ARGS( TEST_CMND_STR "${TEST_CMND_ARRAY}" )
    #PRINT_VAR(TEST_CMND_STR)

    APPEND_STRING_VAR( TEST_SCRIPT_STR
      "\n"
      "SET( TEST_${TEST_CMND_IDX}_CMND ${TEST_CMND_STR} )\n"
      )
    IF (PACKAGE_ADD_ADVANCED_TEST_UNITTEST)
      GLOBAL_SET(PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX}
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
    "SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/utils)\n"
    "\n"
    "INCLUDE(DriveAdvancedTest)\n"
    "\n"
    "DRIVE_ADVANCED_TEST()\n"
    )

  IF (PACKAGE_ADD_ADVANCED_TEST_UNITTEST)
    GLOBAL_SET(PACKAGE_ADD_ADVANCED_TEST_NUM_CMNDS ${NUM_CMNDS})
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(TEST_SCRIPT_STR)
  ENDIF()

  # Write the script file

  SET(TEST_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.cmake")

  IF (ADD_THE_TEST AND NOT PACKAGE_ADD_ADVANCED_TEST_SKIP_SCRIPT)

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("\nWriting file \"${TEST_SCRIPT_FILE}\" ...")
    ENDIF()
  
    FILE( WRITE "${TEST_SCRIPT_FILE}"
      "${TEST_SCRIPT_STR}" )

  ENDIF()

  #
  # F) Set the CTest test to run the new script
  #

  IF (ADD_THE_TEST AND NOT PACKAGE_ADD_ADVANCED_TEST_SKIP_SCRIPT AND NOT PACKAGE_ADD_ADVANCED_TEST_SKIP_SCRIPT_ADD_TEST)

    # Tell CTest to run our script for this test.  Pass the test-time
    # configuration name to the script in the TEST_CONFIG variable.
    ADD_TEST( ${TEST_NAME}
      ${CMAKE_COMMAND} "-DTEST_CONFIG=\${CTEST_CONFIGURATION_TYPE}" -P "${TEST_SCRIPT_FILE}"
      )

    PACKAGE_PRIVATE_ADD_TEST_ADD_LABEL_AND_KEYWORDS(${TEST_NAME})

    #This if clause will set the number of PROCESSORS to reserve during testing
    #to the number requested for the test.
    IF(NUM_PROCS_USED)
      SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES
        PROCESSORS "${NUM_PROCS_USED}")
    ENDIF()
    IF (PARSE_FINAL_PASS_REGULAR_EXPRESSION)
      SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        PASS_REGULAR_EXPRESSION "${PARSE_FINAL_PASS_REGULAR_EXPRESSION}" )
    ELSEIF (PARSE_FINAL_FAIL_REGULAR_EXPRESSION)
      SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        FAIL_REGULAR_EXPRESSION "${PARSE_FINAL_FAIL_REGULAR_EXPRESSION}" )
    ELSE()
      SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
        PASS_REGULAR_EXPRESSION "OVERALL FINAL RESULT: TEST PASSED" )
    ENDIF()

  ENDIF()

ENDFUNCTION()
