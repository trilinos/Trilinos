INCLUDE(PackageAddExecutableAndTest)

IF (WIN32 AND BUILD_SHARED_LIBS)
ELSE()
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    FEVector_test
    SOURCES cxx_main.cpp ExecuteTestProblems.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "tests FAILED"
    DEFINES -DEPETRA_LIB_EXPORTS_MODE
    )
ENDIF()
