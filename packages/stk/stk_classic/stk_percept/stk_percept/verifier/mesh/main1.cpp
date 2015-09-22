#include <iostream>
#include <stk_percept/Percept.hpp>

#include "Verifier.hpp"
#include <stk_percept/PerceptMesh.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_percept/RunEnvironment.hpp>


using namespace stk_classic;
using namespace percept;

#define doMPI 0

int main(int argc,  char **argv)
{
  
#if doMPI
  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) { 
    std::cerr << "MPI_Init FAILED" << std::endl ; 
    std::abort(); 
  } 
#endif

  Verifier vf;
  vf.verify(argc, argv);


#if doMPI
#endif
  //  MPI_Finalize(); 

  return 0;
}
