
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

APPEND_SET(SOURCES
    UnitTestAlgorithmRunner.cpp
    UnitTestCudaMgr.cpp
    UnitTestMain.cpp
    UnitTest_helpers.cpp
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    stk_algsup_unit_tests
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
    )
