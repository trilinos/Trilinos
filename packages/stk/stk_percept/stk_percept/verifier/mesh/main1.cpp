#include <iostream>
#include <stk_percept/Percept.hpp>

#include "Verifier.hpp"
#include <stk_percept/PerceptMesh.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <stk_percept/RunEnvironment.hpp>


using namespace stk;
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

#if 0
  testing::InitGoogleTest(&argc, argv);  
  //  bool result = 0;
  bool result = RUN_ALL_TESTS(); 
  if (!result) return 1;
#endif

  //RunEnvironment::doSierraLoadBalance();

  Verifier vf;
  vf.verify(argc, argv);
  //stk::parallel_machine_finalize();


#if doMPI
#endif
  //  MPI_Finalize(); 

  return 0;
}
