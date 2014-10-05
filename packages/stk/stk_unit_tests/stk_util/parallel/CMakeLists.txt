
INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)
INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_BINARY_DIR}/stk_util)


FILE(GLOB SOURCES *.cpp)
#removing due to dependency on Sierra's unitTestUtils/comparison/stringAndNumberComparisons.h
LIST(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/TimerHowToParallel.cpp)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    util_parallel_UnitTest
    SOURCES ${SOURCES}
    COMM serial mpi
    NUM_MPI_PROCS 1-6
    )
