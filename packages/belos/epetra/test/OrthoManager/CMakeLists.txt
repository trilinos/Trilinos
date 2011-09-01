
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_TSQR)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Epetra_OrthoManager_test
    SOURCES belos_orthomanager_epetra.cpp
    ARGS ""
    COMM serial mpi
    )
ENDIF()
