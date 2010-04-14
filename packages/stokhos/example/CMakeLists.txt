INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES hermite_example.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Example Passed"
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  exp_moment_example
  SOURCES exp_moment_example.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  stieltjes_example
  SOURCES stieltjes_example.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  gram_schmidt_example
  SOURCES gram_schmidt_example.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE(
  recurrence_example
  SOURCES recurrence_example.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  random_field_example
  SOURCES random_field_example.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  sparsity_example
  SOURCES sparsity_example.cpp
  )

IF(Stokhos_ENABLE_TriKota)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    pecos_hermite_example
    SOURCES pecos_hermite_example.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Example Passed"
    NUM_MPI_PROCS 1
  )

  IF(Stokhos_ENABLE_ML)

     PACKAGE_ADD_EXECUTABLE(
       twoD_diffusion_galerkin_example
       SOURCES twoD_diffusion_galerkin_example.cpp
       COMM serial mpi
     )

     PACKAGE_ADD_EXECUTABLE(
       twoD_diffusion_collocation_example
       SOURCES twoD_diffusion_collocation_example.cpp
       COMM serial mpi
     )

  ENDIF()

ENDIF()


