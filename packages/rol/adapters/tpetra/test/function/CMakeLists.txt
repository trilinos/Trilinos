
IF(${PROJECT_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TpetraBoundConstraintInterface
    SOURCES test_01.cpp
    ARGS PrintItAll
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()


