
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/krino/unit_test_utils)

FILE(GLOB SOURCES *.cpp)

# Akri_Unit_Optimize.cpp includes Akri_ROLOptimize.hpp (krino_rol), which sits
# above krino_math_utils.  Under TriBITS the whole krino package is auto-linked,
# so rol is not listed here (doing so is deprecated).  The SierraCMake and bjam
# builds do not auto-link the package, so they list rol explicitly on this test
# target.  krino_lib (external-only) is never used.
TRIBITS_ADD_EXECUTABLE(
    krino_math_utils_utest
    SOURCES ${SOURCES}
    COMM serial mpi
    NOEXEPREFIX NOEXESUFFIX
    TESTONLYLIBS krino_unit_test_utils_lib)

TRIBITS_ADD_TEST(krino_math_utils_utest
    NAME Krino_Math_Utils
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
