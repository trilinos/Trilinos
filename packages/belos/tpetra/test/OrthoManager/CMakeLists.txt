
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tpetra_OrthoManager_test
    SOURCES belos_orthomanager_tpetra.cpp
    ARGS "--ortho=DGKS"
    COMM serial mpi
    )
ENDIF()
