// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <iostream>

#include "MeshTransfer.hpp"

using namespace percept;

int main(int argc,  char **argv)
{
  stk::ParallelMachine comm(stk::parallel_machine_init(&argc, &argv));
  Kokkos::initialize(argc,argv);  

  { // destroy object after run
    MeshTransfer ev(comm);
    ev.run(argc, argv);
  }

  Kokkos::finalize();
  stk::parallel_machine_finalize();

  return 0;
}
