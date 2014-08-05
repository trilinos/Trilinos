
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_BINARY_DIR}/stk_util/stk_util)

FILE(GLOB SOURCES *.cpp)

#removing due to dependence on a header not in Trilinos
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestDistributedIndexWithBulkData.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCreateEdges.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestModificationEnd.cpp")

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    stk_mesh_unit_tests
    SOURCES ${SOURCES}
    DEPLIBS stk_util_use_cases
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    )

