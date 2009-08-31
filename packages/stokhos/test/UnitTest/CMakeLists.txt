
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

IF(Stokhos_ENABLE_Quadrature)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    QuadExpansionUnitTest
    SOURCES Stokhos_QuadExpansionUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    StieltjesUnitTest
    SOURCES Stokhos_StieltjesUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    GramSchmidtUnitTest
    SOURCES Stokhos_GramSchmidtUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Sparse3TensorUnitTest
    SOURCES Stokhos_Sparse3TensorUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ENDIF()
