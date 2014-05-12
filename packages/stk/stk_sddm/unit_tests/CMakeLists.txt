INCLUDE(TribitsAddExecutableAndTest)

#INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_BINARY_DIR}/stk_util)

FILE(GLOB SOURCES *.cpp)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    sddm_unit_tests
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

