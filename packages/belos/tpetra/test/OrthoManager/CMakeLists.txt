
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  IF (Belos_ENABLE_TSQR)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Tpetra_OrthoManager_Tsqr_test
      SOURCES cxx_main_tsqr.cpp 
      ARGS "--ortho=TSQR"
      COMM serial mpi
      )
  ENDIF ()
ENDIF()
