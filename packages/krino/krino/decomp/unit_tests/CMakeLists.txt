
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/krino/unit_test_utils)

FILE(GLOB SOURCES *.cpp)

# Akri_Unit_CDMesh.cpp exercises rebalance_utils::rebalance_mesh / have_parmetis
# and includes Akri_MeshClone.hpp (mesh_clone), which sit above krino_decomp.
# Under TriBITS the whole krino package is auto-linked, so those intra-package
# libs are not listed here (doing so is deprecated).  The SierraCMake and bjam
# builds do not auto-link the package, so they list rebalance_utils/mesh_clone
# explicitly on this test target.  krino_lib (external-only) is never used.
TRIBITS_ADD_EXECUTABLE(
    krino_decomp_utest
    SOURCES ${SOURCES}
    COMM serial mpi
    NOEXEPREFIX NOEXESUFFIX
    TESTONLYLIBS krino_unit_test_utils_lib)

TRIBITS_ADD_TEST(krino_decomp_utest
    NAME Krino_Basic
    ARGS "--gtest_filter=-*Random*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)

TRIBITS_ADD_TEST(krino_decomp_utest
    NAME Krino_Robustness
    ARGS "--gtest_filter=*Random*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
