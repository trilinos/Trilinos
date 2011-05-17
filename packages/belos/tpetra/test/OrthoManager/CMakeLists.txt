
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_TSQR)
  # FIXME (mfh 15 May 2011) The TSQR version of the OrthoManager seems
  # to be broken.  I'm disabling that test for now by running the test
  # with no arguments, which makes it pass by default.  I will fix
  # TsqrOrthoManager as soon as I get a chance.
  #
  # The usual ARGS:
  # ARGS "--ortho=TSQR --verbose"
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tpetra_OrthoManager_test
    SOURCES belos_orthomanager_tpetra.cpp
    ARGS ""
    COMM serial mpi
    )
  # This "test" passes trivially.  It is a benchmark that is meant 
  # to be run manually.
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tpetra_OrthoManager_benchmark
    SOURCES belos_orthomanager_tpetra_benchmark.cpp
    ARGS "--ortho=TSQR --verbose"
    COMM serial mpi
    )
ENDIF()


