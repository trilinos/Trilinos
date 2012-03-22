INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test/src)

IF (WIN32 AND BUILD_SHARED_LIBS)
ELSE()
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    FECrsMatrix_test
    SOURCES cxx_main.cpp ExecuteTestProblems.cpp
    DEPLIBS epetratest
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "tests FAILED"
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    FEmatrix2
    SOURCES fematrix2.cpp
    DEPLIBS epetratest
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "tests FAILED"
    )
ENDIF()
