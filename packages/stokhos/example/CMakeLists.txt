INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE(
  GMRES_example
  SOURCES GMRES.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )
 
TRIBITS_ADD_EXECUTABLE(
  CG_example
  SOURCES CG.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )





TRIBITS_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES hermite_example.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Example Passed"
  NUM_MPI_PROCS 1
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  exp_moment_example
  SOURCES exp_moment_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_example
  SOURCES stieltjes_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_example2
  SOURCES stieltjes_example2.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_example3
  SOURCES stieltjes_example3.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_example4
  SOURCES stieltjes_example4.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_coupled_example
  SOURCES stieltjes_coupled.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  gram_schmidt_example
  SOURCES gram_schmidt_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  gram_schmidt_example2
  SOURCES gram_schmidt_example2.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  recurrence_example
  SOURCES recurrence_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  random_field_example
  SOURCES random_field_example.cpp
  COMM serial mpi
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  sparsity_example
  SOURCES sparsity_example.cpp
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  cijk_nonzeros
  SOURCES cijk_nonzeros.cpp
  DEPLIBS stokhos
  )

TRIBITS_ADD_EXECUTABLE(
  sparsity_slice_example
  SOURCES sparsity_slice_example.cpp
  DEPLIBS stokhos
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(copySparsityScripts
      DEST_FILES   animate_slices.py
      )

TRIBITS_ADD_LIBRARY(
  linear2d_diffusion
  HEADERS twoD_diffusion_ME.hpp
  SOURCES twoD_diffusion_ME.cpp
  TESTONLY
  DEPLIBS stokhos
  )

IF(Stokhos_ENABLE_TriKota)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pecos_hermite_example
    SOURCES pecos_hermite_example.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Example Passed"
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
  )

  IF(Stokhos_ENABLE_AztecOO AND Stokhos_ENABLE_ML)

     TRIBITS_ADD_EXECUTABLE(
       Linear2D_Diffusion_Collocation_Example
       SOURCES linear2d_diffusion_collocation.cpp
       COMM serial mpi
       DEPLIBS stokhos linear2d_diffusion
     )

     IF(Stokhos_ENABLE_Stratimikos)

       TRIBITS_ADD_EXECUTABLE(
         Linear2D_Diffusion_Collocation_Strat_Example
         SOURCES linear2d_diffusion_collocation_strat.cpp
         COMM serial mpi
         DEPLIBS stokhos linear2d_diffusion
       )

       IF(Stokhos_ENABLE_NOX)

         TRIBITS_ADD_EXECUTABLE_AND_TEST(
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

     TRIBITS_ADD_EXECUTABLE_AND_TEST(
       Linear2D_Diffusion_PCE_Example
       SOURCES linear2d_diffusion_pce.cpp
       COMM serial mpi
       DEPLIBS stokhos linear2d_diffusion
       ARGS 1
       NUM_MPI_PROCS 2
       PASS_REGULAR_EXPRESSION "Example Passed!"
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
       Linear2D_Diffusion_PCE_Interlaced_Example
       SOURCES linear2d_diffusion_pce_interlaced.cpp
       COMM serial mpi
       DEPLIBS stokhos linear2d_diffusion
       NUM_MPI_PROCS 2
       PASS_REGULAR_EXPRESSION "Example Passed!"
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(copyInterlacedScripts
      DEST_FILES   check_interlaced_solution.py
      )

      IF(Stokhos_ENABLE_NOX)

        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          Linear2D_Diffusion_PCE_NOX_Example
          SOURCES linear2d_diffusion_pce_nox.cpp
          COMM serial mpi
	  DEPLIBS stokhos linear2d_diffusion
	  ARGS 1
	  NUM_MPI_PROCS 2
	  PASS_REGULAR_EXPRESSION "Example Passed!"
        )

       IF(Stokhos_ENABLE_Stratimikos)
         TRIBITS_ADD_EXECUTABLE(
           Linear2D_Diffusion_PCE_NOX_SG_Solvers
           SOURCES linear2d_diffusion_pce_nox_sg_solvers.cpp
           COMM serial mpi
	   DEPLIBS stokhos linear2d_diffusion
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GMRES_Mean_Based"
	    ARGS "--sg_prec_method=Mean-Based --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

          TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GMRES_AGS"
	    ARGS "--sg_prec_method=Approx-Gauss-Seidel --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_CG_AGS"
	    ARGS "--outer_krylov_method=CG --sg_prec_method=Approx-Gauss-Seidel --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GMRES_GS"
	    ARGS "--sg_prec_method=Gauss-Seidel --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GMRES_AJ"
	    ARGS "--sg_prec_method=Approx-Jacobi --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GMRES_KP"
	    ARGS "--sg_prec_method=Kronecker-Product --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GS"
	    ARGS "--sg_solver=Gauss-Seidel --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_JA"
	    ARGS "--sg_solver=Jacobi --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_LN"
	    ARGS "--rand_field=Log-Normal --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

          TRIBITS_ADD_TEST(
            Linear2D_Diffusion_PCE_NOX_SG_Solvers
            NAME "Linear2D_Diffusion_GSLN"
            ARGS "--sg_solver=Gauss-Seidel --rand_field=Log-Normal --num_spatial_procs=1"
            COMM serial mpi
            NUM_MPI_PROCS 2
            PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GMRES_FA"
	    ARGS "--sg_operator_method=Fully-Assembled --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GMRES_KL"
	    ARGS "--sg_operator_method=KL-Matrix-Free --unnormalize --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )

	  TRIBITS_ADD_TEST(
	    Linear2D_Diffusion_PCE_NOX_SG_Solvers
	    NAME "Linear2D_Diffusion_GMRES_KLR"
	    ARGS "--sg_operator_method=KL-Reduced-Matrix-Free --num_spatial_procs=1"
	    COMM serial mpi
	    NUM_MPI_PROCS 2
	    PASS_REGULAR_EXPRESSION "Example Passed!"
          )
       ENDIF()

        
     ENDIF()  

ENDIF()

IF(Stokhos_ENABLE_Sacado)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    sacado_example
    SOURCES sacado_example.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    sacado_vector_example
    SOURCES sacado_vector_example.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
    )

  IF(Stokhos_ENABLE_KokkosArray)
    
    SET(CUDA_SOURCES "")
    IF (Stokhos_ENABLE_CUDA)
      CUDA_COMPILE(CUDA_SOURCES sacado_mpvector_example_cuda.cu)
    ENDIF()

    TRIBITS_ADD_EXECUTABLE(
      sacado_mpvector_example
      SOURCES sacado_mpvector_example.cpp
              sacado_mpvector_example_host.cpp
	      ${CUDA_SOURCES}
      COMM serial mpi
      DEPLIBS stokhos
      )
  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    division_example
    SOURCES division_example.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Example Passed"
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
    )

  TRIBITS_ADD_EXECUTABLE(
    gram_schmidt_example3
    SOURCES gram_schmidt_example3.cpp
    COMM serial mpi
    DEPLIBS stokhos
    )

  TRIBITS_ADD_EXECUTABLE(
    stieltjes_example5
    SOURCES stieltjes_example5.cpp
    COMM serial mpi
    DEPLIBS stokhos
    )

ENDIF() 

IF(Stokhos_ENABLE_NOX)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    nox_example
    SOURCES nox_example.cpp SimpleME.hpp SimpleME.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS stokhos
    )

ENDIF()

IF (Stokhos_ENABLE_Thrust)

  TRIBITS_ADD_EXECUTABLE(
    cuda_example
    SOURCES cuda_example.cpp
    COMM serial mpi
    DEPLIBS stokhos
  )

  IF(Stokhos_ENABLE_Sacado)

    CUDA_COMPILE(CUDA_SOURCES sacado_cuda_example.cu)
    TRIBITS_ADD_EXECUTABLE(
      sacado_cuda_example
      SOURCES ${CUDA_SOURCES}
      COMM serial mpi
      DEPLIBS stokhos
    )

  ENDIF() 

ENDIF()

IF(Stokhos_ENABLE_Isorropia)

  TRIBITS_ADD_EXECUTABLE(
    cijk_partition
    SOURCES cijk_partition.cpp
    DEPLIBS stokhos
    )

ENDIF()

IF(Stokhos_ENABLE_Experimental)
  ADD_SUBDIRECTORY(experimental)
ENDIF()
