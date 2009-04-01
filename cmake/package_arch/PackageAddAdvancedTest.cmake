
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
#     TEST_0 [EXEC <execTarget0> | CMND <cmndExec0>] ...
#     [TEST_1 [EXEC <execTarget1> | CMND <cmndExec1>] ...]
#     ...
#     [TEST_N [EXEC <execTargetN> | CMND <cmndExecN>] ...]
#     [HOST <host1> <host2> ...]
#     [XHOST <host1> <host2> ...]
#     [COMM [serial] [mpi]]
#     [OVERALL_NUM_PROCS]
#     [FINAL_PASS_REGULAR_EXPRESSION <regex> | FINAL_FAIL_REGULAR_EXPRESSION <regex>]
#     )
#
# Here <testName> is the name of the test (which will have ${PACKAGE_NAME}_
# appended) that will be used to name the output CMake script file as well as
# the CTest test name passed into ADD_TEST(...).
#
# Each and every atomic test or command needs to pass (as defined below) in
# order for the overall test to pass.
#
# Each atomic test case is either a package-built executable or just a basic
# command.  An atomic test or command takes the form:
#
#   TEST_i
#      EXEC <execTarget> [NOEXEPREFIX] [NOEXESUFFIX] | CMND <cmndExec>
#      ARGS <arg1> <arg2> ... <argn>
#      [MESSAGE "<message>"]
#      [WORKING_DIRECTORY <workingDir>]
#      [OUTPUT_FILE <outputFile> [ECHO_OUTPUT_FILE]]
#      [NUM_MPI_PROCS <number>]
#      [PASS_ANY
#        | PASS_REGULAR_EXPRESSION <regex>
#        | PASS_REGULAR_EXPRESSION_ALL "<regex1>;<regex2>;...;<regexn>"
#        | FAIL_REGULAR_EXPRESSION <regex>
#        | STANDARD_PASS_OUTPUT
#        ]
#
# Each test line is either package-built test executable or some general
# command executable:
#
#   EXEC <execTarget>
#
#     If it is a package executable, then you specify just the executable
#     target name as "EXEC <execTarget>".  The value <execTarget> is same
#     string that was passed in as the first argument to
#     PACKAGE_ADD_EXECUTABLE( <execTarget>...) used to define the executable.
#     If this is an MPI build, then <execTarget> will be run with MPI using
#     NUM_MPI_PROCS <number> or OVERALL_NUM_MPI_PROCS <number> (if
#     NUM_MPI_PROCS is not set for this test case..  If the number of maximum
#     MPI processes allowed is less than this number of MPI processes, then
#     the test will *not* be run.  If NOEXEPREFIX is specified, then
#     ${PACKAGE_NAME}_ will not be added the the beginning.  If NOEXESUFFIX is
#     specified, then '.exe' will not be added to the end.
#
#   CMND <cmndExec>
#
#     If it is a general command, then it takes the form "CMND <cmndExec>".
#     When this is a general command, then the command is run as is without
#     MPI.  In this case, MPI will not be used to run the executable even when
#     configured in MPI mode (i.e. TPL_ENABLE_MPI=ON).  Note that EXEC
#     <execTarget> is basically equivalent to CMND <cmndExec> when NOEXEPREFIX
#     and NOEXESUFFIX are specified.  In this case, you can pass in
#     <execTarget> to any command you would like and it will get run with MPI
#     in MPI mode just link any other command.
#
# Other miscellaneous arguments:
#
#   OUTPUT_FILE <outputFile>
#
#     If specified, <outputFile> is the name of a file that standard out will
#     be captured into (relative to the current working directory).  Capturing
#     the output to a file allows a later command to read the file and perform
#     some type of test (e.g. like a diff or something).
#
#   ECHO_OUTPUT_FILE
#
#     If specified, then the contents of the file <outputFile> given in the
#     OUTPUT_FILE argument will be echoed to the screen after the command has
#     finished.
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
#     from this command in some way.  Therefore, you would typically use
#     "OUTPUT_FILE <outputFile>" when using this.
#
#   PASS_REGULAR_EXPRESSION <regex>
#
#     If specified, the test command 'i' will be assumed to pass if it matches
#     the given regular expression.  Otherwise, it is assumed to fail.
#
#   FAIL_REGULAR_EXPRESSION <regex>
#
#     If specified, the test command 'i' will be assumed to fail if it matches
#     the given regular expression.  Otherwise, it is assumed to pass.
#
#   STANDARD_PASS_OUTPUT
#
#     If specified, the test command 'i' will be assumed to pass if the string
#     expression "Final Result: PASSED" is found in the ouptut for the test.
#
# ???
#

FUNCTION(PACKAGE_ADD_ADVANCED_TEST TEST_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_ADD_ADVANCED_TEST: ${TEST_NAME}\n")
  ENDIF()

  SET(TEST_NAME ${PACKAGE_NAME}_${TEST_NAME_IN})

  #
  # A) Parse the overall arguments and figure out how many tests
  # comands we will have
  #

  # Allow for a maximum of 10 (0 through 9) test commands
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
     "${TEST_IDX_LIST};COMM"
     #options
     ""
     ${ARGN}
     )

  #
  # E) Determine if we will add the serial or MPI tests based on input COMM
  # and TPL_ENABLE_MPI
  #

  PACKAGE_PROCESS_COMM_ARGS(ADD_SERIAL_TEST  ADD_MPI_TEST  ${PARSE_COMM})
  IF (NOT ADD_SERIAL_TEST AND NOT ADD_MPI_TEST)
    RETURN()
  ENDIF()

  #
  # B) Build the test script
  #

  SET(TEST_SCRIPT_STR "")

  APPEND_STRING_VAR( TEST_SCRIPT_STR
    "\n"
    "#\n"
    "# This is a CMake script and must be run as \"cmake -P SCRIPT_NAME\"\n"
    "#\n"
    "\n"
    "#\n"
    "# Variables\n"
    "#\n"
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
       "EXEC;CMND;ARGS;OUTPUT_FILE;NUM_MPI_PROCS"
       #options
       "NOEXEPREFIX;NOEXESUFFIX;PASS_ANY;STANDARD_PASS_OUTPUT"
       ${PARSE_TEST_${TEST_CMND_IDX}}
       )
    # Write the command

    IF (PARSE_EXEC)
      LIST( LENGTH PARSE_EXEC PARSE_EXEC_LEN )
      IF (NOT PARSE_EXEC_LEN EQUAL 1)
        MESSAGE(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} EXEC = '${PARSE_EXEC}'"
          " must be a single name.  To add arguments use ARGS <arg1> <arg2> ...." )
      ENDIF()
      PACAKGE_ADD_TEST_GET_EXE_BINARY_NAME( "${PARSE_EXEC}"
        ${PARSE_NOEXEPREFIX} ${PARSE_NOEXESUFFIX} EXE_BINARY_NAME )
      SET( TEST_CMND_ARRAY "./${EXE_BINARY_NAME}" )
    ELSEIF (PARSE_CMND)
      LIST( LENGTH PARSE_CMND PARSE_CMND_LEN )
      IF (NOT PARSE_CMND_LEN EQUAL 1)
        MESSAGE(SEND_ERROR "Error, TEST_${TEST_CMND_IDX} CMND = '${PARSE_CMND}'"
          " must be a single command.  To add arguments use ARGS <arg1> <arg2> ...." )
      ENDIF()
      SET( TEST_CMND_ARRAY ${PARSE_CMND} )
    ELSE()
      MESSAGE( FATAL_ERROR
        "Must have at least EXEC or CMND for TEST_${TEST_CMND_IDX}" )
    ENDIF()

    IF (PARSE_ARGS)
      JOIN_EXEC_PROCESS_SET_ARGS( ARGS_STR ${PARSE_ARGS} )
      SET( TEST_CMND_ARRAY "${TEST_CMND_ARRAY} ${ARGS_STR}" )
    ENDIF()

    APPEND_STRING_VAR( TEST_SCRIPT_STR
      "\n"
      "SET( TEST_CMND_${TEST_CMND_IDX} ${TEST_CMND_ARRAY} )\n"
      )
    IF (PACKAGE_ADD_ADVANCED_TEST_UNITTEST)
      GLOBAL_SET(PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX}
        "${TEST_CMND_ARRAY}" )
    ENDIF()

    IF (PARSE_OUTPUT_FILE)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_OUTPUT_FILE_${TEST_CMND_IDX} \"${PARSE_OUTPUT_FILE}\" )\n"
        )
    ENDIF()

  ENDFOREACH()

  # ToDo: Verify that TEST_${MAX_NUM_TEST_CMND_IDX}+1 does *not* exist!

  APPEND_STRING_VAR( TEST_SCRIPT_STR
    "\n"
    "SET(NUM_CMNDS ${NUM_CMNDS})\n"
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

  # Write script the file

  SET(TEST_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.cmake")

  IF (NOT PACKAGE_ADD_ADVANCED_TEST_SKIP_SCRIPT)

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("\nWriting file \"${TEST_SCRIPT_FILE}\" ...")
    ENDIF()
  
    FILE( WRITE "${TEST_SCRIPT_FILE}"
      "${TEST_SCRIPT_STR}" )

  ENDIF()

  #
  # C) Set the CTest test to run the new script
  #

  IF (NOT PACKAGE_ADD_ADVANCED_TEST_SKIP_SCRIPT)

    ADD_TEST( ${TEST_NAME}
      ${CMAKE_COMMAND} -P "${TEST_SCRIPT_FILE}"
      )

    SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
      PASS_REGULAR_EXPRESSION "FINAL RESULT: TEST PASSED" )

  ENDIF()

ENDFUNCTION()
