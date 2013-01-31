INCLUDE(TribitsAddExecutableAndTest)

IF(Stokhos_ENABLE_Tpetra AND Stokhos_ENABLE_Belos AND Stokhos_ENABLE_Ifpack2 AND Stokhos_ENABLE_MueLu)

     TRIBITS_ADD_EXECUTABLE(
       Linear2D_Diffusion_Example
       SOURCES linear2d_diffusion.cpp 
               twoD_diffusion_problem_tpetra.hpp 
	       twoD_diffusion_problem_tpetra_def.hpp
       COMM serial mpi
       DEPLIBS stokhos
      )

    TRIBITS_ADD_EXECUTABLE(
       Linear2D_Diffusion_IPCE_Example
       SOURCES linear2d_diffusion_pce.cpp 
               twoD_diffusion_problem_tpetra.hpp 
	       twoD_diffusion_problem_tpetra_def.hpp
       COMM serial mpi
       DEPLIBS stokhos muelu xpetra
      )
    TRIBITS_ADD_EXECUTABLE(
       Linear2D_Diffusion_IPCE_Ifpack2_Example
       SOURCES linear2d_diffusion_pce_ifpack2.cpp 
               twoD_diffusion_problem_tpetra.hpp 
	       twoD_diffusion_problem_tpetra_def.hpp
       COMM serial mpi
       DEPLIBS stokhos
      )
     TRIBITS_ADD_TEST(
            Linear2D_Diffusion_IPCE_Ifpack2_Example
            NAME "Linear2D_Diffusion_IPCE_Ifpack2_Example_Direct"
            ARGS "--num_mesh=8 --division_method=Direct"
            COMM serial mpi
            NUM_MPI_PROCS 1 
            PASS_REGULAR_EXPRESSION "Example Passed!"
          )
       TRIBITS_ADD_TEST(
            Linear2D_Diffusion_IPCE_Ifpack2_Example
            NAME "Linear2D_Diffusion_IPCE_Ifpack2_Example_CG_Diag"
            ARGS "--num_mesh=8 --division_method=CG --divprec_method=Diag"
            COMM serial mpi
            NUM_MPI_PROCS 1
            PASS_REGULAR_EXPRESSION "Example Passed!"
          )
      TRIBITS_ADD_TEST(
            Linear2D_Diffusion_IPCE_Ifpack2_Example
            NAME "Linear2D_Diffusion_IPCE_Ifpack2_Example_CG_Schur"
            ARGS "--num_mesh=8 --division_method=CG --divprec_method=Schur"
            COMM serial mpi
            NUM_MPI_PROCS 1
            PASS_REGULAR_EXPRESSION "Example Passed!"
          )




ENDIF()

