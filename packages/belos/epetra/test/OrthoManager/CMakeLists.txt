
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_TSQR)
  # FIXME (mfh 05 Nov 2010) TSQR is currently the only OrthoManager
  # subclass that passes this test (using the Euclidean inner
  # product).  This is because none of DGKS, IMGS, or ICGS currently
  # do the right thing for rank-deficient inputs.  Thus, we only build
  # this test if TSQR is enabled.
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Epetra_OrthoManager_test
    SOURCES belos_orthomanager_epetra.cpp
    ARGS "--ortho=TSQR"
    COMM serial mpi
    )
ENDIF()
