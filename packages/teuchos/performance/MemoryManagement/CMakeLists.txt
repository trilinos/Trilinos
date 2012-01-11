
# Required modules
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)

#
# RCP_PerforamnceTests
#

TRIBITS_ADD_EXECUTABLE(
  RCP_PerformanceTests
  COMM serial mpi  
  SOURCES
    RCP_Performance_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )

# This form of the test just makes sure that it runs without crashing
# basically
TRIBITS_ADD_TEST(
  RCP_PerformanceTests
  NUM_MPI_PROCS 1
  COMM serial mpi
  POSTFIX_AND_ARGS_0 base
    --show-test-details=ALL --max-array-size=100000
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS)

  TRIBITS_ADD_TEST(
    RCP_PerformanceTests
    NUM_MPI_PROCS 1
    COMM serial mpi
    POSTFIX_AND_ARGS_0 performance
      --show-test-details=ALL --max-array-size=100000 --rel-cpu-speed=1e+3
      --max-rcp-create-destroy-ratio=1.05
      --max-rcp-raw-adjust-ref-count-ratio=20.0 --max-rcp-sp-adjust-ref-count-ratio=1.7
      --max-rcp-raw-obj-access-ratio=1.02
    HOST godel.sandia.gov gabriel.sandia.gov
    STANDARD_PASS_OUTPUT
    )

ENDIF()

#
# Array_PerforamnceTests
#

TRIBITS_ADD_EXECUTABLE(
  Array_PerformanceTests
  COMM serial mpi  
  SOURCES
    Array_Performance_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )

# This form of the test just makes sure that it runs without crashing
# basically
TRIBITS_ADD_TEST(
  Array_PerformanceTests
  NUM_MPI_PROCS 1
  COMM serial mpi
  POSTFIX_AND_ARGS_0 base
    --show-test-details=ALL --max-array-size=100000
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS)

  TRIBITS_ADD_TEST(
    Array_PerformanceTests
    NUM_MPI_PROCS 1
    COMM serial mpi
    POSTFIX_AND_ARGS_0 performance
      --details=ALL --max-array-size=100000 --rel-cpu-speed=1e+3
      --max-array-bracket-ratio=1.02 --max-array-iter-ratio=1.02
      --max-arrayrcp-self-iter-ratio=1.5
    HOST godel.sandia.gov gabriel.sandia.gov
    STANDARD_PASS_OUTPUT
    )

ENDIF()
