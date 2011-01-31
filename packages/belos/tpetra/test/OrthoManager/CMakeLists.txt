
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_TSQR)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tpetra_OrthoManager_test
    SOURCES belos_orthomanager_tpetra.cpp
    ARGS "--ortho=TSQR --verbose"
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


