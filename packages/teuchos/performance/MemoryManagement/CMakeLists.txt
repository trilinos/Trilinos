
# Required modules
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)

#
# RCP_PerforamnceTests
#

PACKAGE_ADD_EXECUTABLE(
  RCP_PerformanceTests
  COMM serial mpi  
  SOURCES
    RCP_Performance_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  )

# This form of the test just makes sure that it runs without crashing
# basically
PACKAGE_ADD_TEST(
  RCP_PerformanceTests
  NAME RCP_PerformanceTests_base
  NUM_MPI_PROCS 1
  COMM serial mpi
    ARGS
      "--show-test-details=ALL --max-array-size=16384"
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS)

  PACKAGE_ADD_TEST(
    RCP_PerformanceTests
    NAME RCP_PerformanceTests_performance
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS
      "--show-test-details=ALL --max-array-size=262144 --rel-cpu-speed=1e+2 --max-rcp-create-destroy-ratio=1.50 --max-rcp-adjust-ref-count-ratio=1.10"
    HOST godel.sandia.gov gabriel.sandia.gov
    STANDARD_PASS_OUTPUT
    )

ENDIF()

#
# Array_PerforamnceTests
#

PACKAGE_ADD_EXECUTABLE(
  Array_PerformanceTests
  COMM serial mpi  
  SOURCES
    Array_Performance_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  )

# This form of the test just makes sure that it runs without crashing
# basically
PACKAGE_ADD_TEST(
  Array_PerformanceTests
  NAME Array_PerformanceTests_base
  NUM_MPI_PROCS 1
  COMM serial mpi
    ARGS
      "--show-test-details=ALL --max-array-size=16384"
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS)

  PACKAGE_ADD_TEST(
    Array_PerformanceTests
    NAME Array_PerformanceTests_performance
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS
      "--details=ALL --max-array-size=262144 --rel-cpu-speed=1e+3 --max-array-bracket-ratio=1.02 --max-arrayrcp-self-iter-ratio=1.5"
    HOST godel.sandia.gov gabriel.sandia.gov
    STANDARD_PASS_OUTPUT
    )

ENDIF()
