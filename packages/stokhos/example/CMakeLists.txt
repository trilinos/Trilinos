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

     PACKAGE_ADD_EXECUTABLE(
       twoD_diffusion_collocation_example_rcg-v2
       SOURCES twoD_diffusion_collocation_example-rcg-v2.cpp
       COMM serial mpi
     )

     PACKAGE_ADD_EXECUTABLE(
       twoD_diffusion_collocation_log-normal_example
       SOURCES twoD_diffusion_collocation_log-normal.cpp
       COMM serial mpi
     )

     PACKAGE_ADD_EXECUTABLE(
       Linear2D_Diffusion_PCE_Example
       SOURCES linear2d_diffusion_pce.cpp
               twoD_diffusion_ME.hpp twoD_diffusion_ME.cpp
       COMM serial mpi
      )

     PACKAGE_ADD_EXECUTABLE(
       Linear2D_Diffusion_PCE_NOX_Example
       SOURCES linear2d_diffusion_pce_nox.cpp
               twoD_diffusion_ME.hpp twoD_diffusion_ME.cpp
       COMM serial mpi
      )

     PACKAGE_ADD_EXECUTABLE(
       Linear2D_Diffusion_PCE_NOX_SG_Solvers
       SOURCES linear2d_diffusion_pce_nox_sg_solvers.cpp
               twoD_diffusion_ME.hpp twoD_diffusion_ME.cpp
       COMM serial mpi
      )

     PACKAGE_ADD_EXECUTABLE(
       Linear2D_Diffusion_PCE_NOX_Belos
       SOURCES linear2d_diffusion_pce_nox_belos.cpp
               twoD_diffusion_ME.hpp twoD_diffusion_ME.cpp
       COMM serial mpi
      )  

  ENDIF()

ENDIF()

IF(Stokhos_ENABLE_Sacado)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    sacado_example
    SOURCES sacado_example.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

ENDIF() 

IF(Stokhos_ENABLE_NOX)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    nox_example
    SOURCES nox_example.cpp SimpleME.hpp SimpleME.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

ENDIF()
