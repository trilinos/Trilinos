INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES hermite_example.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Example Passed"
  NUM_MPI_PROCS 1
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE(
  exp_moment_example
  SOURCES exp_moment_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE(
  stieltjes_example
  SOURCES stieltjes_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE(
  stieltjes_example2
  SOURCES stieltjes_example2.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE(
  stieltjes_coupled_example
  SOURCES stieltjes_coupled.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

PACKAGE_ADD_EXECUTABLE(
  gram_schmidt_example
  SOURCES gram_schmidt_example.cpp
  COMM serial mpi
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

PACKAGE_ADD_EXECUTABLE(
  sparsity_slice_example
  SOURCES sparsity_slice_example.cpp
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

       IF(Stokhos_ENABLE_NOX)

         PACKAGE_ADD_EXECUTABLE_AND_TEST(
           Linear2D_Diffusion_PCE_MPNI_Example
           SOURCES linear2d_diffusion_pce_mpni.cpp
           COMM serial mpi
	   DEPLIBS stokhos linear2d_diffusion
	   ARGS 1
	   NUM_MPI_PROCS 2
	   PASS_REGULAR_EXPRESSION "Example Passed!"
         )

       ENDIF()

     ENDIF()
  ENDIF()
ENDIF()

IF(Stokhos_ENABLE_AztecOO AND Stokhos_ENABLE_ML)

     PACKAGE_ADD_EXECUTABLE_AND_TEST(
       Linear2D_Diffusion_PCE_Example
       SOURCES linear2d_diffusion_pce.cpp
       COMM serial mpi
       DEPLIBS stokhos linear2d_diffusion
       ARGS 1
       NUM_MPI_PROCS 2
       PASS_REGULAR_EXPRESSION "Example Passed!"
      )

      IF(Stokhos_ENABLE_NOX)

        PACKAGE_ADD_EXECUTABLE_AND_TEST(
          Linear2D_Diffusion_PCE_NOX_Example
          SOURCES linear2d_diffusion_pce_nox.cpp
          COMM serial mpi
	  DEPLIBS stokhos linear2d_diffusion
	  ARGS 1
	  NUM_MPI_PROCS 2
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
	  ARGS "--sg_prec_method=Mean-Based --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

        PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_AGS"
	  ARGS "--sg_prec_method=Approx-Gauss-Seidel --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_CG_AGS"
	  ARGS "--outer_krylov_method=CG --sg_prec_method=Approx-Gauss-Seidel --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_GS"
	  ARGS "--sg_prec_method=Gauss-Seidel --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_AJ"
	  ARGS "--sg_prec_method=Approx-Jacobi --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_KP"
	  ARGS "--sg_prec_method=Kronecker-Product --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GS"
	  ARGS "--sg_solver=Gauss-Seidel --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GS"
	  ARGS "--sg_solver=Jacobi --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_LN"
	  ARGS "--rand_field=Log-Normal --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_FA"
	  ARGS "--sg_operator_method=Fully-Assembled --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_KL"
	  ARGS "--sg_operator_method=KL-Matrix-Free --unnormalize --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

	PACKAGE_ADD_TEST(
	  Linear2D_Diffusion_PCE_NOX_SG_Solvers
	  NAME "Linear2D_Diffusion_GMRES_KLR"
	  ARGS "--sg_operator_method=KL-Reduced-Matrix-Free --num_spatial_procs=1"
	  COMM serial mpi
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
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
    DEPLIBS stokhos
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    sacado_vector_example
    SOURCES sacado_vector_example.cpp
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


