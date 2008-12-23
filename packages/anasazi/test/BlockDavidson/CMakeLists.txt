
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)

IF(Teuchos_ENABLE_COMPLEX)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    BlockDavidson_complex_lap_test
    SOURCES cxx_main_complex_lap.cpp
    ARGS --verbose
    COMM serial mpi
    )
ENDIF()
