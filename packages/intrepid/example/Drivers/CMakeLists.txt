INCLUDE(PackageAddExecutableAndTest)

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Example_01
    SOURCES example_01.cpp
    DEPLIBS epetra epetraext
    ARGS "10 10 10 false 1.0 10.0 0.0 1.0 -1.0 1.0 -1.0 1.0"
    COMM serial
#    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

ENDIF()
