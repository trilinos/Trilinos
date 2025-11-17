

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_OrthoManager_test
  SOURCES belos_orthomanager_tpetra.cpp
  ARGS "--ortho=ICGS --verbose"
  ARGS "--ortho=IMGS --verbose"
  ARGS "--ortho=DGKS --verbose"
  COMM serial mpi
  )
# This "test" passes trivially with the given arguments.  That is
# because it is a benchmark that is meant to be run manually, with
# carefully chosen arguments.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_OrthoManager_benchmark
  SOURCES belos_orthomanager_tpetra_benchmark.cpp
  ARGS "--ortho=ICGS --verbose"
  COMM serial mpi
  )

IF (${PACKAGE_NAME}_ENABLE_TSQR)
  TRIBITS_ADD_TEST(
    Tpetra_OrthoManager_test
    NAME Tpetra_OrthoManager_test_TSQR
    ARGS "--ortho=TSQR --verbose"
    COMM serial mpi
    )
  # This "test" passes trivially with the given arguments.  That is
  # because it is a benchmark that is meant to be run manually, with
  # carefully chosen arguments.
  TRIBITS_ADD_TEST(
    Tpetra_OrthoManager_benchmark
    NAME Tpetra_OrthoManager_benchmark_TSQR
    ARGS "--ortho=TSQR --verbose"
    COMM serial mpi
    )
ENDIF()
