
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

IF(Stokhos_ENABLE_Quadrature)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    UnitTest
    SOURCES Stokhos_UnitTest.cpp Stokhos_StieltjesUnitTest.cpp Stokhos_QuadExpansionUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ENDIF()
