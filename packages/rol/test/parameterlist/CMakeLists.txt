#GET_PROPERTY( PTR_STRING GLOBAL PROPERTY PTR_IMPL )
IF( ROL_ParameterList STREQUAL "simple")
  message("Enabling tests of ROL's native ParameterList implementation")
  ROL_ADD_EXECUTABLE_AND_TEST(
    ParameterList 
    SOURCES test_01.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()
