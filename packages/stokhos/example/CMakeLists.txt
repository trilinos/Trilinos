INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES hermite_example.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Example Passed"
  NUM_MPI_PROCS 1
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  exp_moment_example
  SOURCES exp_moment_example.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  stieltjes_example
  SOURCES stieltjes_example.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  gram_schmidt_example
  SOURCES gram_schmidt_example.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE(
  recurrence_example
  SOURCES recurrence_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE(
  random_field_example
  SOURCES random_field_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE(
  sparsity_example
  SOURCES sparsity_example.cpp
  DEPLIBS stokhos
  )

PACKAGE_ADD_LIBRARY(
  linear2d_diffusion
  HEADERS twoD_diffusion_ME.hpp
  SOURCES twoD_diffusion_ME.cpp
  TESTONLY
  DEPLIBS stokhos
  )

IF(Stokhos_ENABLE_TriKota)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    pecos_hermite_example
    SOURCES pecos_hermite_example.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Example Passed"
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
  )

  IF(Stokhos_ENABLE_AztecOO AND Stokhos_ENABLE_ML)

     PACKAGE_ADD_EXECUTABLE(
       Linear2D_Diffusion_Collocation_Example
       SOURCES linear2d_diffusion_collocation.cpp
       COMM serial mpi
       DEPLIBS stokhos linear2d_diffusion
     )

     IF(Stokhos_ENABLE_Stratimikos)

       PACKAGE_ADD_EXECUTABLE(
         Linear2D_Diffusion_Collocation_Strat_Example
         SOURCES linear2d_diffusion_collocation_strat.cpp
         COMM serial mpi
         DEPLIBS stokhos linear2d_diffusion
       )

     ENDIF()

     PACKAGE_ADD_EXECUTABLE_AND_TEST(
       Linear2D_Diffusion_PCE_Example
       SOURCES linear2d_diffusion_pce.cpp
       COMM serial mpi
       DEPLIBS stokhos linear2d_diffusion
       NUM_MPI_PROCS 4
       PASS_REGULAR_EXPRESSION "Example Passed!"
      )

      IF(Stokhos_ENABLE_NOX)

        PACKAGE_ADD_EXECUTABLE_AND_TEST(
          Linear2D_Diffusion_PCE_NOX_Example
          SOURCES linear2d_diffusion_pce_nox.cpp
          COMM serial mpi
	  DEPLIBS stokhos linear2d_diffusion
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

       PACKAGE_ADD_EXECUTABLE(
         Linear2D_Diffusion_PCE_NOX_SG_Solvers
         SOURCES linear2d_diffusion_pce_nox_sg_solvers.cpp
         COMM serial mpi
	 DEPLIBS stokhos linear2d_diffusion
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_Mean_Based"
	  ARGS "--sg_prec_method=Mean-Based"
	  COMM serial mpi
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

        PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_AGS"
	  ARGS "--sg_prec_method=Approx-Gauss-Seidel"
	  COMM serial mpi
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_GS"
	  ARGS "--sg_prec_method=Gauss-Seidel"
	  COMM serial mpi
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_AJ"
	  ARGS "--sg_prec_method=Approx-Jacobi"
	  COMM serial mpi
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_KP"
	  ARGS "--sg_prec_method=Kronecker-Product"
	  COMM serial mpi
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GS"
	  ARGS "--sg_solver=Gauss-Seidel"
	  COMM serial mpi
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GS"
	  ARGS "--sg_solver=Jacobi"
	  COMM serial mpi
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_LN"
	  ARGS "--rand_field=Log-Normal"
	  COMM serial mpi
	  NUM_MPI_PROCS 4
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )
	  

     ENDIF()  

  ENDIF()

ENDIF()

IF(Stokhos_ENABLE_Sacado)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    sacado_example
    SOURCES sacado_example.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
    )

ENDIF() 

IF(Stokhos_ENABLE_NOX)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    nox_example
    SOURCES nox_example.cpp SimpleME.hpp SimpleME.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
    )

ENDIF()

IF (Stokhos_ENABLE_Thrust)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    cuda_example
    SOURCES cuda_example.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
  )

  IF(Stokhos_ENABLE_Sacado)

    CUDA_COMPILE(CUDA_SOURCES sacado_cuda_example.cu)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      sacado_cuda_example
      SOURCES ${CUDA_SOURCES}
      ARGS -v
      COMM serial mpi
      NUM_MPI_PROCS 1
      DEPLIBS stokhos
    )

  ENDIF() 

ENDIF()


