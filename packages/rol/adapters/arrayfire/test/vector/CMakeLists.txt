
IF(ROL_ENABLE_ArrayFireCPU)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ArrayFireVectorInterface
    SOURCES test_01.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()


