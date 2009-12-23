
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  LegendreBasisUnitTest
  SOURCES Stokhos_LegendreBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  HermiteBasisUnitTest
  SOURCES Stokhos_HermiteBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DerivExpansionUnitTest
  SOURCES Stokhos_DerivExpansionUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  QuadExpansionUnitTest
  SOURCES Stokhos_QuadExpansionUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

IF(TPL_ENABLE_ForUQTK)
 
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ForUQTKExpansionUnitTest
    SOURCES Stokhos_ForUQTKExpansionUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ENDIF()

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

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ExponentialRandomFieldUnitTest
  SOURCES Stokhos_ExponentialRandomFieldUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
