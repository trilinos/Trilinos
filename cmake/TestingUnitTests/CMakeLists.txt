
INCLUDE(MultilineSet)
INCLUDE(GenerateScriptFileTarget)


FUNCTION(CREATE_TESTING_UNITTEST_TEST_CASE TEST_NAME)

  MULTILINE_SET(TEST_CMND
    "#!/bin/bash\n"
    "if [ ! -e ${TEST_NAME} ]\n"
    "  then mkdir ${TEST_NAME}\n"
    "fi\n"
    "cd ${TEST_NAME}\n"
    "if [ -e CMakeCache.txt ]\n"
    "  then rm CMakeCache.txt\n"
    "fi\n"
    "\"${CMAKE_COMMAND}\""
    "   -DTrilinos_INVOKE_TESTING_UNIT_TESTS:BOOL=ON"
    "  ${TRILINOS_HOME_DIR}\n"
    )
  
  SET(TEST_SCRIPT_FILE "${TEST_NAME}.sh")

  GENERATE_SCRIPT_FILE_TARGET(${TEST_SCRIPT_FILE} ${TEST_CMND}
    ${TRILINOS_HOME_DIR}/CMakeLists.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    )
  
  ADD_TEST( ${TEST_NAME} ${TEST_SCRIPT_FILE} )
  
  SET_TESTS_PROPERTIES( ${TEST_NAME}
    PROPERTIES PASS_REGULAR_EXPRESSION "Final UnitTests Result: PASSED" )

ENDFUNCTION()


CREATE_TESTING_UNITTEST_TEST_CASE(TrilinosUnitTests_TestingUnitTest)
