

ASSERT_DEFINED(Belos_ENABLE_Triutils)

IF (Belos_ENABLE_Triutils AND Belos_ENABLE_TSQR)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Epetra_OrthoManager_test
    SOURCES belos_orthomanager_epetra.cpp
    ARGS ""
    COMM serial mpi
    )
ENDIF()
