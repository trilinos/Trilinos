IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_BlockDavidson_test
    SOURCES cxx_main.cpp
    ARGS --verbose
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()
