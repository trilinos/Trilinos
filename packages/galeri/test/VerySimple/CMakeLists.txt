
IF (${PACKAGE_NAME}_ENABLE_Epetra)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    VerySimple
    SOURCES cxx_main.cpp
    ARGS -v
    COMM serial
    )

ENDIF()
