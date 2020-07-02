

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Epetra_OrthoManager_test
    SOURCES belos_orthomanager_epetra.cpp
    ARGS "--ortho=ICGS --verbose --debug"
         "--ortho=DGKS --verbose --debug"
         "--ortho=IMGS --verbose --debug"
    COMM serial mpi
    )

ENDIF()
