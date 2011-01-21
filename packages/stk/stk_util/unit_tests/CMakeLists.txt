
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

ADD_SUBDIRECTORY(parallel)

FILE(GLOB SOURCES *.cpp diag/*.cpp environment/*.cpp util/*.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    util_UnitTest
    SOURCES ${SOURCES}
    DEPLIBS stk_util_parallel stk_util_diag
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )


