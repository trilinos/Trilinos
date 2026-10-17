
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/krino/unit_test_utils)

FILE(GLOB SOURCES *.cpp)

TRIBITS_ADD_EXECUTABLE(
    krino_refinement_utest
    SOURCES ${SOURCES}
    COMM serial mpi
    NOEXEPREFIX NOEXESUFFIX
    TESTONLYLIBS krino_unit_test_utils_lib)

TRIBITS_ADD_TEST(krino_refinement_utest
    NAME Krino_Refinement
    ARGS "--gtest_filter=-*fuzz*:*Refinement.performance*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)

TRIBITS_ADD_TEST(krino_refinement_utest
    NAME Krino_Fuzz
    ARGS "--gtest_filter=*fuzz*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)

TRIBITS_ADD_TEST(krino_refinement_utest
    NAME Krino_Refinement_Performance
    ARGS "--gtest_filter=*Refinement.performance*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
