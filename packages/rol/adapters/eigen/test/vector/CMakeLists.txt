
IF(ROL_ENABLE_Eigen)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    EigenVectorInterface
    SOURCES test_01.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()


