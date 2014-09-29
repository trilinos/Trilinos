
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_BINARY_DIR}/stk_util/stk_util)

FILE(GLOB SOURCES *.cpp)
#removing due to dependence on geometry toolkit which is not available to Trilinos
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCoarseSearch.cpp")

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    search_unit_tests
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

