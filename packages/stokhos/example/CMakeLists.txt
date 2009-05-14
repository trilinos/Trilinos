INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES hermite_example.cpp
  COMM serial mpi
  )

IF(Stokhos_ENABLE_Quadrature)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    product_basis_example
    SOURCES product_basis_example.cpp
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    exp_example
    SOURCES exp_example.cpp
    COMM serial mpi
    )

ENDIF()

PACKAGE_ADD_EXECUTABLE(
  dakota_pce
  SOURCES dakota_pce.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  dakota_exp
  SOURCES dakota_pce.cpp
  COMM serial mpi
  )
