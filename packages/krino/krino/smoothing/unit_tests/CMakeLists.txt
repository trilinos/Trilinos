
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/krino/unit_test_utils)

FILE(GLOB SOURCES *.cpp)

# The smoothing tests include Akri_ROLOptimize.hpp (krino_rol); the smoothing lib
# itself does not depend on rol.  Under TriBITS the whole krino package is
# auto-linked, so rol is not listed here (doing so is deprecated).  The
# SierraCMake and bjam builds do not auto-link the package, so they list rol
# explicitly on this test target.  krino_lib (external-only) is never used.
TRIBITS_ADD_EXECUTABLE(
    krino_smoothing_utest
    SOURCES ${SOURCES}
    COMM serial mpi
    NOEXEPREFIX NOEXESUFFIX
    TESTONLYLIBS krino_unit_test_utils_lib)

TRIBITS_ADD_TEST(krino_smoothing_utest
    NAME Krino_Smoothing
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
