
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/krino/unit_test_utils)

FILE(GLOB SOURCES *.cpp)

# Akri_Unit_FastMarching.cpp includes Akri_PostProcess.hpp (post_process), which
# sits above krino_eikonal.  Under TriBITS the whole krino package is auto-linked,
# so post_process is not listed here (doing so is deprecated).  The SierraCMake
# and bjam builds do not auto-link the package, so they list post_process
# explicitly on this test target.  krino_lib (external-only) is never used.
TRIBITS_ADD_EXECUTABLE(
    krino_eikonal_utest
    SOURCES ${SOURCES}
    COMM serial mpi
    NOEXEPREFIX NOEXESUFFIX
    TESTONLYLIBS krino_unit_test_utils_lib)

TRIBITS_ADD_TEST(krino_eikonal_utest
    NAME Krino_Eikonal
    ARGS "--gtest_filter=-HamiltonJacobi.*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)

TRIBITS_ADD_TEST(krino_eikonal_utest
    NAME Krino_HamiltonJacobi
    ARGS "--gtest_filter=HamiltonJacobi.*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
