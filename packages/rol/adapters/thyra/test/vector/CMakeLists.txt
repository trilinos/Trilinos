
IF(${PROJECT_NAME}_ENABLE_Thyra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ThyraInterface
    SOURCES test_01.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    ) 
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ScaledThyraVector
    SOURCES test_02.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    ) 
ENDIF()
