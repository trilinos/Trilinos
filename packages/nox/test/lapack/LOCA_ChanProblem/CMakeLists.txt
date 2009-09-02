
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK AND NOX_ENABLE_LOCA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_ChanContinuation
    SOURCES ChanContinuation.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS noxtestutils localapacktestproblems
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_ChanTPContinuation
    SOURCES ChanTPContinuation.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS noxtestutils localapacktestproblems
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_ChanHomotopy
    SOURCES ChanHomotopy.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS noxtestutils localapacktestproblems
    PASS_REGULAR_EXPRESSION "All tests passed"
    )
  
ENDIF()
