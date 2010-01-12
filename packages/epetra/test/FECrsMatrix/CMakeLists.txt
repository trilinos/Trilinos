INCLUDE(PackageAddExecutableAndTest)

IF (WIN32 AND BUILD_SHARED_LIBS)
ELSE()
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    FECrsMatrix_test
    SOURCES cxx_main.cpp ExecuteTestProblems.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "tests FAILED"
    )
ENDIF()
