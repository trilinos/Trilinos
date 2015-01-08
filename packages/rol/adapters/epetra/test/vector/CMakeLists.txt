INCLUDE(TribitsAddExecutableAndTest)

IF(${PROJECT_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    EpetraMultiVectorInterface
    SOURCES test_01.cpp
    DEPLIBS epetra
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()

