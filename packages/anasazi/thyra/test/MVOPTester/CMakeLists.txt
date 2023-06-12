

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MVOPTesterThyra_test
    SOURCES cxx_main.cpp
    ARGS -v
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()
