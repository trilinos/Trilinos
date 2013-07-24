INCLUDE(TribitsAddExecutableAndTest)

IF(Stokhos_ENABLE_Tpetra AND Stokhos_ENABLE_Belos AND Stokhos_ENABLE_Ifpack2
AND KokkosClassic_ENABLE_Thrust AND KokkosClassic_ENABLE_ThreadPool)

     TRIBITS_ADD_EXECUTABLE(
       Linear2D_Diffusion_Example
       SOURCES linear2d_diffusion.cpp 
               twoD_diffusion_problem_tpetra.hpp 
	       twoD_diffusion_problem_tpetra_def.hpp
       COMM serial mpi
       DEPLIBS stokhos tpetra
      )

ENDIF()

