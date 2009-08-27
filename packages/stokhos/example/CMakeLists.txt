INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES hermite_example.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  recurrence_example
  SOURCES recurrence_example.cpp
  COMM serial mpi
  )

IF(TPL_ENABLE_Dakota)

  IF(Stokhos_ENABLE_ML)
     PACKAGE_ADD_EXECUTABLE(
       twoD_diffusion_example
       SOURCES twoD_diffusion_example.cpp
       COMM serial mpi
     )

     PACKAGE_ADD_EXECUTABLE(
       sparse_grid_example
       SOURCES twoD_diffusion_sparse_grid.cpp
       COMM serial mpi
     )
  ENDIF()

ENDIF()

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

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    moment_example
    SOURCES moment_example.cpp
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    stieltjes_example
    SOURCES stieltjes_example.cpp
    COMM serial mpi
    )

   PACKAGE_ADD_EXECUTABLE_AND_TEST(
    stieltjes_dependent_example
    SOURCES stieltjes_dependent_example.cpp
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    gram_schmidt_example
    SOURCES gram_schmidt_example.cpp
    COMM serial mpi
    )

ENDIF()

PACKAGE_ADD_EXECUTABLE(
  dakota_pce
  SOURCES dakota_pce.cpp
  COMM serial mpi
  )
