// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <percept/Percept.hpp>

#include "Verifier.hpp"
#include <percept/PerceptMesh.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Kokkos_Core.hpp>

#include <percept/RunEnvironment.hpp>


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
  Kokkos::initialize(argc, argv);

  Verifier vf;
  vf.verify(argc, argv);


  Kokkos::finalize();
#if doMPI
#endif
  //  MPI_Finalize(); 

  return 0;
}
