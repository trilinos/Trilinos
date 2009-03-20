
INCLUDE(ParseVariableArguments)
INCLUDE(AppendStringVar)
INCLUDE(Join)
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
#     <TestName>
#     TEST_0 [ EXEC <ExecTarget> | CMND <cmndExec>] ...
#     TEST_1 [ EXEC <ExecTarget> | CMND <cmndExec>] ...
#     ...
#     TEST_N [ EXEC <ExecTarget> | CMND <cmndExec>] ...
#     [ COMM [serial] [mpi] ]
#     [ OVERALL_NUM_MPI_PROCS <number> ]
#     [ HOST ]
#     [ XHOST ]
#     )
#
# Here 'TestName' is the name of the test (which will have ${PACKAGE_NAME}_
# appended) that will be used to name the output CMake script file as well as
# the CTest test name passed into ADD_TEST(...).
#
# Each and every atomic test or command needs to pass (as defined below) in
# order for the overall test to pass.
#
# Each atomic test case is either a package-built executable or just a basic
# command.  An atomic test or command takes the form:
#
#   TEST_i [EXEC <ExecTarget> | CMND <cmndExec>]
#      ARGS <arg1> <arg2> ... <argn>
#      [ NOEXEPREFIX ]
#      [ NOEXESUFFIX ]
#      [ OUTPUT_FILE <outputFile> ]
#      [ NUM_MPI_PROCS <number> ]
#      [ PASS_ANY
#          | PASS_REGULAR_EXPRESSION <regex>
#          | FAIL_REGULAR_EXPRESSION <regex>
#          | STANDARD_PASS_OUTPUT ]
#
# Each test test line is either package-built test executable or some general
# command executable.  If it is a package executable, then you specify just
# the executable target name as "EXEC ExecTarget".  This is the same string
# that was passed in as the first argument to the PACKAGE_ADD_EXECUTABLE(...) 
# function.  If it is a general command, then it takes the form "CMND
# <cmndExec>".  If this is an MPI build, then ExecTarget will be run with MPI
# using NUM_MPI_PROCS <number>.  If the number of maximum MPI processes
# allowed is less than this number of MPI processes, then the test will *not*
# be run.  When this is a general command, then the command is run as is
# without MPI.
#

# If OUTPUT_FILE is given, the next argment <outputFile> is the name of a file
# that standard out will be captured into (relative to the current working
# directory).  Capturing the output to a file allows a later command to read
# the file and perform some type of test (e.g. like a diff or something).


#
# By default, an atomic test line is assumed to pass if the executable returns
# a non-zero value.  However, a test case can also be defined to pass based
# on:
#
#   PASS_ANY
#      If specified, the test command 'i' will be assumed to pass reguardless
#      of the return value or any other output.  This would be used when a
#      command that is to follow will determine pass or fail based on output
#      from this command in some way.
#
#   PASS_REGULAR_EXPRESSION <regex>
#
#     If specified, the test command 'i' will be assumed to pass if it matches
#     the given regular expression.  Otherwise, it is assumed to fail.
#
#   FAIL_REGULAR_EXPRESSION <regex>
#
#     The test command 'i' will be assumed to fail if it matches the given
#     regular expression.  Otherwise, it is assumed to pass.
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

  SET(MAX_NUM_TEST_CMND_IDX 9)

  SET(TEST_IDX_LIST "")
  FOREACH( TEST_CMND_IDX RANGE ${MAX_NUM_TEST_CMND_IDX})
    LIST( APPEND TEST_IDX_LIST TEST_${TEST_CMND_IDX} )
  ENDFOREACH()
  #PRINT_VAR(TEST_IDX_LIST)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "${TEST_IDX_LIST}"
     #options
     ""
     ${ARGN}
     )

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
      MESSAGE( FATAL_ERROR "EXEC not implmeneted yet!" )
    ELSEIF (PARSE_CMND)
      SET( TEST_CMND_ARRAY ${PARSE_CMND} )
    ELSE()
      MESSAGE( FATAL_ERROR
        "Must have at least EXEC or CMND for TEST_${TEST_CMND_IDX}" )
    ENDIF()

    IF (PARSE_ARGS)
      JOIN( ARGS_STR " " ${PARSE_ARGS} )
      LIST( APPEND TEST_CMND_ARRAY "${ARGS_STR}" )
    ENDIF()

    APPEND_STRING_VAR( TEST_SCRIPT_STR
      "\n"
      "SET( TEST_CMND_${TEST_CMND_IDX} ${TEST_CMND_ARRAY} )\n"
      )

    IF (PARSE_OUTPUT_FILE)
      APPEND_STRING_VAR( TEST_SCRIPT_STR
        "\n"
        "SET( TEST_OUTPUT_FILE_${TEST_CMND_IDX} \"${PARSE_OUTPUT_FILE}\" )\n"
        )
    ENDIF()

  ENDFOREACH()

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

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(TEST_SCRIPT_STR)
  ENDIF()

  # Write script the file

  SET(TEST_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.cmake")

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nWriting file \"${TEST_SCRIPT_FILE}\" ...")
  ENDIF()

  FILE( WRITE "${TEST_SCRIPT_FILE}"
    "${TEST_SCRIPT_STR}" )

  #
  # C) Set the CTest test to run the new script
  #

  ADD_TEST( ${TEST_NAME}
    ${CMAKE_COMMAND} -P "${TEST_SCRIPT_FILE}"
    )

  SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
    PASS_REGULAR_EXPRESSION "FINAL RESULT: TEST PASSED" )

ENDFUNCTION()
