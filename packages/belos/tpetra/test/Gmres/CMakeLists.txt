
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils)
  # This "test" passes trivially.  It is a benchmark that is meant 
  # to be run manually.
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tpetra_Gmres_test
    SOURCES belos_gmres_tpetra.cpp
    ARGS "--verbose"
    COMM serial mpi
    )
ENDIF()


