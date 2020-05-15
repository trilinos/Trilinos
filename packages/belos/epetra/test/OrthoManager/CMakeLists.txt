

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Epetra_OrthoManager_test
    SOURCES belos_orthomanager_epetra.cpp
    ARGS "--ortho=ICGS --verbose"
         "--ortho=DGKS --verbose"
         "--ortho=IMGS --verbose"
    COMM serial mpi
    )

ENDIF()
