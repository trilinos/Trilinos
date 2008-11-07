
# Required modules
INCLUDE(Trilinos_Add_Executable)
INCLUDE(Trilinos_Add_Test)

TRILINOS_ADD_EXECUTABLE(
  MemoryManagementPerformanceTests
  COMM serial mpi  
  SOURCES
    RCP_Performance_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  )

# This form of the test just makes sure that it runs without crashing
# basically
TRILINOS_ADD_TEST(
  MemoryManagementPerformanceTests
  NAME MemoryManagementPerformanceTests_base
  NUM_MPI_PROCS 1
  COMM serial mpi
    ARGS
      "--show-test-details=ALL --max-array-size=16384"
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS)

  TRILINOS_ADD_TEST(
    MemoryManagementPerformanceTests
    NAME MemoryManagementPerformanceTests_performance
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS
      "--show-test-details=ALL --max-array-size=16384 --rel-cpu-speed=1e+2 --max-rcp-create-destroy-ratio=1.40"
    HOST godel.sandia.gov
    STANDARD_PASS_OUTPUT
    )

ENDIF()
