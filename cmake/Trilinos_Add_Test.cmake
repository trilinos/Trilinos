

# 2008/07/09: rabartl: ToDo:
#
# (*) Support optional DIRECTORY argument
#
# (*) Support multiple ARGS keywords
#
# (*) Support an optional POSTFIX argument for naming tests with different
# ARGS keywords

#INCLUDE(Trilinos_Add_Executable)
INCLUDE(Parse_Variable_Arguments)

FUNCTION (TRILINOS_ADD_TEST EXE_NAME)
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     PARSE   #prefix
     "NUM_MPI_PROCS;DIRECTORY;KEYWORDS;COMM;ARGS;NAME;PASS_REGULAR_EXPRESSION;HOST;XHOST;FAIL_REGULAR_EXPRESSION"  #lists
     "DESEND_INTO_DIR"   #options
     ${ARGN} )

  IF(PARSE_ARGS)
    LIST(LENGTH PARSE_ARGS NUM_PARSE_ARGS)
    #MESSAGE(STATUS "NUM_PARSE_ARGS = ${NUM_PARSE_ARGS}")
  ENDIF()

  IF("${VERBOSE_CONFIGURE}" STREQUAL  "ON")
    MESSAGE("")
    MESSAGE("TRILINOS_ADD_TEST: EXE_NAME = ${EXE_NAME}")
  ENDIF()
  
  #
  # B) Do not add any tests if host or xhost says not to
  #
  
  IF(NOT PARSE_XHOST)
    SET (PARSE_XHOST NONE)
  ENDIF()    
  LIST (FIND PARSE_XHOST ${TRILINOS_HOSTNAME} INDEX_OF_HOSTNAME_IN_XHOST_LIST)           
  IF (NOT ${INDEX_OF_HOSTNAME_IN_XHOST_LIST} EQUAL -1)
    RETURN()
  ENDIF()

  IF(NOT PARSE_HOST)
    SET (PARSE_HOST ${TRILINOS_HOSTNAME})
  ENDIF()  
  LIST (FIND PARSE_HOST ${TRILINOS_HOSTNAME} INDEX_OF_HOSTNAME_IN_HOST_LIST)                 
  IF (${INDEX_OF_HOSTNAME_IN_HOST_LIST} EQUAL -1)
    RETURN()
  ENDIF()

  #
  # C) Set the name and path of the binary that will be run
  #
  
  SET(EXE_BINARY_NAME "${EXE_NAME}${CMAKE_EXECUTABLE_SUFFIX}")
  #MESSAGE("TRILINOS_ADD_TEST: ${EXE_NAME}: EXE_BINARY_NAME = ${EXE_BINARY_NAME}")

  IF(PARSE_NAME)
    SET(TEST_NAME "${PROJECT_NAME}_${PARSE_NAME}")
  ELSE()
    SET(TEST_NAME "${PROJECT_NAME}_${EXE_NAME}")  
  ENDIF()

  IF(PARSE_DIRECTORY)
    SET(EXECUTABLE_PATH "./${PARSE_DIRECTORY}/${EXE_BINARY_NAME}")
  ELSE()
    SET(EXECUTABLE_PATH "./${EXE_BINARY_NAME}")
  ENDIF()
  #MESSAGE("TRILINOS_ADD_TEST: ${EXE_NAME}: EXECUTABLE_PATH = ${EXECUTABLE_PATH}")

  #
  # D) Append keywords to the name of the test
  #
  
  IF(PARSE_KEYWORDS)
    FOREACH(KEYWORD ${PARSE_KEYWORDS})
      SET(TEST_NAME ${TEST_NAME}_${KEYWORD})
    ENDFOREACH()
  ENDIF()
  
  SET(ADDED_THE_TEST OFF)
    
  IF(TRILINOS_ENABLE_MPI)
    SET(NP)
    SET(NUM_PROCS_USED 1)
    IF(PARSE_NUM_MPI_PROCS)
      IF(${PARSE_NUM_MPI_PROCS} MATCHES [0-9]+-[0-9]+)
        STRING(REGEX REPLACE "([0-9]+)-([0-9]+)" "\\1" MIN_NP ${PARSE_NUM_MPI_PROCS} )
        STRING(REGEX REPLACE "([0-9]+)-([0-9]+)" "\\2" MAX_NP ${PARSE_NUM_MPI_PROCS} )
        IF(${MIN_NP} LESS ${MPIEXEC_MAX_NUMPROCS} AND  ${MAX_NP} GREATER ${MPIEXEC_MAX_NUMPROCS} )
          SET(NUM_PROCS_USED ${MPIEXEC_MAX_NUMPROCS})
        ELSEIF(${MIN_NP} EQUAL ${MPIEXEC_MAX_NUMPROCS})
          SET(NUM_PROCS_USED ${MIN_NP})
        ELSEIF(${MAX_NP} EQUAL ${MPIEXEC_MAX_NUMPROCS})
          SET(NUM_PROCS_USED ${MAX_NP})
        ELSEIF(${MAX_NP} LESS ${MPIEXEC_MAX_NUMPROCS})
          SET(NUM_PROCS_USED ${MAX_NP})
        ELSE()
          # The number of available processor is outside the given range
          # so the test should not be run.
          RETURN()
        ENDIF()
      ELSE()
        IF(${PARSE_NUM_MPI_PROCS} LESS ${MPIEXEC_MAX_NUMPROCS})
          SET(NUM_PROCS_USED ${PARSE_NUM_MPI_PROCS})
        ELSE()
          SET(NUM_PROCS_USED ${MPIEXEC_MAX_NUMPROCS})
        ENDIF()
      ENDIF()
    ENDIF()

    SET(NP ${MPI_NUMPROCS_FLAG} ${NUM_PROCS_USED})
    
    #
    # E) Add the tests
    #
   
    IF (NOT PARSE_COMM)
      # If no COMM is given assume we will add the test
      SET(DO_MPI_INDEX 0)
    ELSE()
      # Else, if COMM is defined we have to find 'mpi'
      LIST (FIND PARSE_COMM "mpi" DO_MPI_INDEX)
    ENDIF()

    IF(NOT ${DO_MPI_INDEX} EQUAL -1)

      SET(TEST_NAME "${TEST_NAME}_MPI_${NUM_PROCS_USED}")
      
      IF(PARSE_ARGS)

        SET(COUNTER 0)

        FOREACH(PARSE_ARG ${PARSE_ARGS})

          IF(${NUM_PARSE_ARGS} EQUAL 1)
            SET(TEST_NAME_COUNTER "${TEST_NAME}")
          ELSE()
            SET(TEST_NAME_COUNTER "${TEST_NAME}_${COUNTER}")
          ENDIF()
          IF("${VERBOSE_CONFIGURE}" STREQUAL "ON")
            MESSAGE(STATUS "TEST_NAME = ${TEST_NAME_COUNTER}")
          ENDIF()
          
          #This is a little bit of a hack
          #If the argument string has multiple arguments then the white space will need 
          #to replaced by a semicolin.  If this is not done the add_test command will
          #add a slash to each white space in the argument string.
          STRING(REPLACE " " ";" MYARG ${PARSE_ARG}) 
          ADD_TEST(${TEST_NAME_COUNTER} ${MPI_EXECUTABLE} ${NP} ${EXECUTABLE_PATH} ${MYARG})
          SET(ADDED_THE_TEST ON)    
          
          IF (PARSE_PASS_REGULAR_EXPRESSION)
            SET_TESTS_PROPERTIES(${TEST_NAME_COUNTER} PROPERTIES PASS_REGULAR_EXPRESSION
              ${PARSE_PASS_REGULAR_EXPRESSION})
          ENDIF()
  
          IF (PARSE_FAIL_REGULAR_EXPRESSION)
            SET_TESTS_PROPERTIES(${TEST_NAME_COUNTER} PROPERTIES FAIL_REGULAR_EXPRESSION
            ${PARSE_FAIL_REGULAR_EXPRESSION})
          ENDIF()

          MATH(EXPR COUNTER ${COUNTER}+1 )

        ENDFOREACH()

      ELSE()

        IF("${VERBOSE_CONFIGURE}" STREQUAL  "ON")
          MESSAGE(STATUS "TEST_NAME = ${TEST_NAME}")
        ENDIF()

        ADD_TEST(${TEST_NAME} ${MPI_EXECUTABLE} ${NP} ${EXECUTABLE_PATH} )
        SET(ADDED_THE_TEST ON)    
        
      ENDIF()

    ENDIF()

  ELSE()
    
    IF (NOT PARSE_COMM)
      # If no COMM is given assume we will add the test
      SET(DO_SERIAL_INDEX 0)
    ELSE()
      # Else, if COMM is defined we have to find 'serial'
      LIST (FIND PARSE_COMM "serial" DO_SERIAL_INDEX)
    ENDIF()
    
    IF(NOT ${DO_SERIAL_INDEX} EQUAL -1)

      IF(PARSE_ARGS)

        SET(COUNTER 0)

        FOREACH(PARSE_ARG ${PARSE_ARGS})

          IF(${NUM_PARSE_ARGS} EQUAL 1)
            SET(TEST_NAME_COUNTER "${TEST_NAME}")
          ELSE()
            SET(TEST_NAME_COUNTER "${TEST_NAME}_${COUNTER}")
          ENDIF()
          IF("${VERBOSE_CONFIGURE}" STREQUAL  "ON")
            MESSAGE(STATUS "TEST_NAME = ${TEST_NAME_COUNTER}")
          ENDIF()
       
          # See above about this hack
          STRING(REPLACE " " ";" MYARG ${PARSE_ARG})
          ADD_TEST(${TEST_NAME_COUNTER} ${EXECUTABLE_PATH} ${MYARG})
          SET(ADDED_THE_TEST ON)    
          
          IF (PARSE_PASS_REGULAR_EXPRESSION)
            SET_TESTS_PROPERTIES( ${TEST_NAME_COUNTER}
              PROPERTIES PASS_REGULAR_EXPRESSION ${PARSE_PASS_REGULAR_EXPRESSION})
          ENDIF()
  
          IF (PARSE_FAIL_REGULAR_EXPRESSION)
            SET_TESTS_PROPERTIES( ${TEST_NAME_COUNTER}
              PROPERTIES FAIL_REGULAR_EXPRESSION ${PARSE_FAIL_REGULAR_EXPRESSION})
          ENDIF()

          MATH(EXPR COUNTER ${COUNTER}+1 )
        
        ENDFOREACH()

      ELSE()

        IF("${VERBOSE_CONFIGURE}" STREQUAL  "ON")
           MESSAGE(STATUS "TEST_NAME = ${TEST_NAME}")
        ENDIF()

        ADD_TEST(${TEST_NAME} ${EXECUTABLE_PATH} )
        SET(ADDED_THE_TEST ON)    

      ENDIF()

    ENDIF()

  ENDIF()

  # 2008/07/09: rabartl: ToDo: Above, create a macho called
  # ???ITEM_EXITS_IN_LIST??(...) to simplify logic!
    
  IF (PARSE_PASS_REGULAR_EXPRESSION AND ${ADDED_THE_TEST} AND NOT PARSE_ARGS)
   SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION
     ${PARSE_PASS_REGULAR_EXPRESSION})
  ENDIF()
  
  IF (PARSE_FAIL_REGULAR_EXPRESSION AND ADDED_THE_TEST AND NOT PARSE_ARGS)
   SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION
     ${PARSE_FAIL_REGULAR_EXPRESSION})
  ENDIF()

  
ENDFUNCTION(TRILINOS_ADD_TEST)

