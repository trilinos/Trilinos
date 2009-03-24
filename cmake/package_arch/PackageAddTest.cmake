INCLUDE(PackageAddTestHelpers)


#
# Add a test or a set of tests for a single executable.
#
# PACKAGE_ADD_TEST(
#   <execName>
#   [ NOEXEPREFIX ]
#   [ NOEXESUFFIX ]
#   [ NAME <testName> ]
#   [ DIRECTORY <directory> ]
#   [ ARGS "<arg1> <arg2> ..." "<arg3> <arg4> ..." ... ]
#   [ COMM [serial] [mpi] ]
#   [ NUM_MPI_PROCS <numProcs> ]
#   [ HOST <host1> <host2> ... ]
#   [ XHOST <host1> <host2> ... ]
#   [ HOSTTYPE <hosttype1> <hosttype2> ... ]
#   [ XHOSTTYPE <hosttype1> <hosttype2> ... ]
#   [ STANDARD_PASS_OUTPUT
#     | PASS_REGULAR_EXPRESSION "<regex1>;<regex2>;..." 
#     | FAIL_REGULAR_EXPRESSION "<regex1>;<regex2>;..." ]
#   )
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
#     created by the function PACKAGE_ADD_EXECUTABLE(<execName> ...).
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
#   DIRECTORY <directory>
#
#     If specified, then the executable is assumed to be in the directory
#     given by relative <directory>.
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

FUNCTION(PACKAGE_ADD_TEST EXE_NAME)
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "DIRECTORY;KEYWORDS;COMM;NUM_MPI_PROCS;ARGS;NAME;HOST;XHOST;HOSTTYPE;XHOSTTYPE;PASS_REGULAR_EXPRESSION;FAIL_REGULAR_EXPRESSION"
     #options
     "NOEXEPREFIX;NOEXESUFFIX;STANDARD_PASS_OUTPUT"
     ${ARGN}
     )

  IF (PARSE_ARGS)
    LIST(LENGTH PARSE_ARGS NUM_PARSE_ARGS)
  ELSE()
    SET(PARSE_ARGS " ")
    SET(NUM_PARSE_ARGS 1)
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("PACKAGE_ADD_TEST: EXE_NAME = ${EXE_NAME}")
  ENDIF()
  
  #
  # B) Do not add any tests if host/hostname or xhost/xhostname says not to
  #

  IF (NOT PARSE_HOST)
    SET (PARSE_HOST ${${PROJECT_NAME}_HOSTNAME})
  ENDIF()  
  LIST (FIND PARSE_HOST ${${PROJECT_NAME}_HOSTNAME} INDEX_OF_HOSTNAME_IN_HOST_LIST)                 
  IF (${INDEX_OF_HOSTNAME_IN_HOST_LIST} EQUAL -1)
    RETURN()
  ENDIF()
  
  IF (NOT PARSE_XHOST)
    SET (PARSE_XHOST NONE)
  ENDIF()    
  LIST (FIND PARSE_XHOST ${${PROJECT_NAME}_HOSTNAME} INDEX_OF_HOSTNAME_IN_XHOST_LIST)           
  IF (NOT ${INDEX_OF_HOSTNAME_IN_XHOST_LIST} EQUAL -1)
    RETURN()
  ENDIF()

  IF (NOT PARSE_HOSTTYPE)
    SET(PARSE_HOSTTYPE ${CMAKE_HOST_SYSTEM_NAME})
  ENDIF()
  LIST (FIND PARSE_HOSTTYPE ${CMAKE_HOST_SYSTEM_NAME} INDEX_OF_HOSTSYSTEMNAME_IN_HOSTTYPE_LIST)
  IF (${INDEX_OF_HOSTSYSTEMNAME_IN_HOSTTYPE_LIST} EQUAL -1)
    RETURN()
  ENDIF()

  IF (NOT PARSE_XHOSTTYPE)
    SET(PARSE_XHOSTTYPE NONE)
  ENDIF()
  LIST (FIND PARSE_XHOSTTYPE ${CMAKE_HOST_SYSTEM_NAME} INDEX_OF_HOSTSYSTEMNAME_IN_XHOSTTYPE_LIST)
  IF (NOT ${INDEX_OF_HOSTSYSTEMNAME_IN_XHOSTTYPE_LIST} EQUAL -1)
    RETURN()
  ENDIF()

  #
  # C) Set the name and path of the binary that will be run
  #
  
  SET(EXE_BINARY_NAME "${EXE_NAME}")
  IF (NOT PARSE_NOEXESUFFIX)
    SET(EXE_BINARY_NAME "${EXE_NAME}${CMAKE_EXECUTABLE_SUFFIX}")
  ENDIF()
  IF(DEFINED PACKAGE_NAME AND NOT PARSE_NOEXEPREFIX)
    SET(EXE_BINARY_NAME ${PACKAGE_NAME}_${EXE_BINARY_NAME})
  ENDIF()

  #MESSAGE("PACKAGE_ADD_TEST: ${EXE_NAME}: EXE_BINARY_NAME = ${EXE_BINARY_NAME}")

  IF(PARSE_NAME)
    SET(TEST_NAME "${PACKAGE_NAME}_${PARSE_NAME}")
  ELSE()
    SET(TEST_NAME "${PACKAGE_NAME}_${EXE_NAME}")  
  ENDIF()

  IF(PARSE_DIRECTORY)
    SET(EXECUTABLE_PATH "${PARSE_DIRECTORY}/${EXE_BINARY_NAME}")
  ELSE()
    SET(EXECUTABLE_PATH "${EXE_BINARY_NAME}")
  ENDIF()

  IF (NOT IS_ABSOLUTE ${EXECUTABLE_PATH})
    SET(EXECUTABLE_PATH "./${EXECUTABLE_PATH}")
  ENDIF()

  #MESSAGE("PACKAGE_ADD_TEST: ${EXE_NAME}: EXECUTABLE_PATH = ${EXECUTABLE_PATH}")

  #
  # D) Append keywords to the name of the test
  #
  
  IF(PARSE_KEYWORDS)
    FOREACH(KEYWORD ${PARSE_KEYWORDS})
      SET(TEST_NAME ${TEST_NAME}_${KEYWORD})
    ENDFOREACH()
  ENDIF()

  #
  # E) Determine if we will add the serial or MPI tests based on input COMM
  # and TPL_ENABLE_MPI
  #

  IF (PARSE_COMM)
    SET(ADD_MPI_TEST OFF)
    SET(ADD_SERIAL_TEST OFF)
    FOREACH(COMM ${PARSE_COMM})
      IF (COMM STREQUAL "mpi")
        SET(ADD_MPI_TEST ON)
      ELSEIF(COMM STREQUAL "serial")
        SET(ADD_SERIAL_TEST ON)
      ELSE()
        MESSAGE(SEND_ERROR "Error, the COMM value '${COMM}' is not valid!.  Only 'mpi' and 'serial' are allowed.")
      ENDIF()
    ENDFOREACH()
  ELSE()
    SET(ADD_MPI_TEST ON)
    SET(ADD_SERIAL_TEST ON)
  ENDIF()

  IF (TPL_ENABLE_MPI)
    SET(ADD_SERIAL_TEST OFF)
  ELSE()
    SET(ADD_MPI_TEST OFF)
  ENDIF()

  #
  # F) Get teh MPI options
  #
    
  IF(TPL_ENABLE_MPI)

    SET(NUM_PROCS_USED ${MPI_EXEC_MAX_NUMPROCS})
    IF(PARSE_NUM_MPI_PROCS)
      IF(${PARSE_NUM_MPI_PROCS} MATCHES [0-9]+-[0-9]+)
        STRING(REGEX REPLACE "([0-9]+)-([0-9]+)" "\\1" MIN_NP ${PARSE_NUM_MPI_PROCS} )
        STRING(REGEX REPLACE "([0-9]+)-([0-9]+)" "\\2" MAX_NP ${PARSE_NUM_MPI_PROCS} )
        IF(${MIN_NP} LESS ${MPI_EXEC_MAX_NUMPROCS} AND  ${MAX_NP} GREATER ${MPI_EXEC_MAX_NUMPROCS} )
          SET(NUM_PROCS_USED ${MPI_EXEC_MAX_NUMPROCS})
        ELSEIF(${MIN_NP} EQUAL ${MPI_EXEC_MAX_NUMPROCS})
          SET(NUM_PROCS_USED ${MIN_NP})
        ELSEIF(${MAX_NP} EQUAL ${MPI_EXEC_MAX_NUMPROCS})
          SET(NUM_PROCS_USED ${MAX_NP})
        ELSEIF(${MAX_NP} LESS ${MPI_EXEC_MAX_NUMPROCS})
          SET(NUM_PROCS_USED ${MAX_NP})
        ELSE()
          # The number of available processor is outside the given range
          # so the test should not be run.
          RETURN()
        ENDIF()
      ELSEIF(${PARSE_NUM_MPI_PROCS} MATCHES [0-9]+,[0-9]+)
        MESSAGE(SEND_ERROR "The test ${TEST_NAME} can not be added yet"
          " because it we do not yet support the form of"
          " NUM_MPI_PROCS=${PARSE_NUM_MPI_PROCS}") 
      ELSE()
        IF(${PARSE_NUM_MPI_PROCS} LESS ${MPI_EXEC_MAX_NUMPROCS})
          SET(NUM_PROCS_USED ${PARSE_NUM_MPI_PROCS})
        ELSE()
          SET(NUM_PROCS_USED ${MPI_EXEC_MAX_NUMPROCS})
        ENDIF()
      ENDIF()
    ENDIF()

  ENDIF()
    
  #
  # G) Add the tests
  #

  IF(ADD_MPI_TEST)

    SET(TEST_NAME "${TEST_NAME}_MPI_${NUM_PROCS_USED}")
    
    SET(COUNTER 0)

    FOREACH(PARSE_ARG ${PARSE_ARGS})

      IF(${NUM_PARSE_ARGS} EQUAL 1)
        SET(TEST_NAME_COUNTER "${TEST_NAME}")
      ELSE()
        SET(TEST_NAME_COUNTER "${TEST_NAME}_${COUNTER}")
      ENDIF()
      IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "TEST_NAME = ${TEST_NAME_COUNTER}")
      ENDIF()

      CONVERT_CMND_ARG_STRING_TO_ADD_TEST_ARG_ARRAY(${PARSE_ARG} INARGS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(INARGS)
      ENDIF()

      PACKAGE_ADD_TEST_ADD_TEST( ${TEST_NAME_COUNTER}
        ${MPI_EXEC}
        ${MPI_EXEC_PRE_NUMPROCS_FLAGS}
        ${MPI_EXEC_NUMPROCS_FLAG} ${NUM_PROCS_USED}
        ${MPI_EXEC_POST_NUMPROCS_FLAGS}
        ${EXECUTABLE_PATH}
        ${INARGS}
        )
      
      PACKAGE_PRIVATE_ADD_TEST_SET_PASS_PROPERTY(${TEST_NAME_COUNTER})

      SET_PROPERTY(TEST ${TEST_NAME_COUNTER} APPEND PROPERTY
        LABELS ${PACKAGE_NAME})

      IF(PARSE_KEYWORDS)
        SET_PROPERTY(TEST ${TEST_NAME_COUNTER} APPEND PROPERTY
          LABELS ${PARSE_KEYWORDS})
      ENDIF()

      MATH(EXPR COUNTER ${COUNTER}+1 )

    ENDFOREACH()

  ENDIF()

  IF(ADD_SERIAL_TEST)

    SET(COUNTER 0)

    FOREACH(PARSE_ARG ${PARSE_ARGS})

      IF(${NUM_PARSE_ARGS} EQUAL 1)
        SET(TEST_NAME_COUNTER "${TEST_NAME}")
      ELSE()
        SET(TEST_NAME_COUNTER "${TEST_NAME}_${COUNTER}")
      ENDIF()
      IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "TEST_NAME = ${TEST_NAME_COUNTER}")
      ENDIF()
    
      CONVERT_CMND_ARG_STRING_TO_ADD_TEST_ARG_ARRAY(${PARSE_ARG} INARGS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(INARGS)
      ENDIF()

      PACKAGE_ADD_TEST_ADD_TEST( ${TEST_NAME_COUNTER} ${EXECUTABLE_PATH} ${INARGS})

      PACKAGE_PRIVATE_ADD_TEST_SET_PASS_PROPERTY(${TEST_NAME_COUNTER})

      IF (NOT PACKAGE_ADD_TEST_ADD_TEST_SKIP)
        SET_PROPERTY(TEST ${TEST_NAME_COUNTER} APPEND PROPERTY
          LABELS ${PACKAGE_NAME})
        IF(PARSE_KEYWORDS)
          SET_PROPERTY(TEST ${TEST_NAME_COUNTER} APPEND PROPERTY
            LABELS ${PARSE_KEYWORDS})
        ENDIF()
      ENDIF()

      MATH(EXPR COUNTER ${COUNTER}+1 )
    
    ENDFOREACH()

  ENDIF()
  
ENDFUNCTION()

