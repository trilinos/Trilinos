INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_BINARY_DIR}/stk_util/stk_util)

FILE(GLOB SOURCES *.cpp)
#disabling due to circular dependency.
LIST(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestSearchWithPeriodicBC.cpp)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    search_util_unit_tests
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

