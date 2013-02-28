
INCLUDE(TribitsAddExecutableAndTest)

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    GeneralizedDavidsonThyra_test
    SOURCES cxx_main.cpp
    ARGS --verbose
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()
