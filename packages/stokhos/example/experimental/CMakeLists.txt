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
               BelosPCETpetraAdapter.hpp
	       tpetra_pce_instantiations.cpp
	       ifpack2_pce_instantiations.cpp
       COMM serial mpi
       DEPLIBS stokhos muelu muelu-xpetra
      )

ENDIF()

