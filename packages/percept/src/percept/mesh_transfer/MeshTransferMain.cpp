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
  MeshTransfer ev;
  ev.run(argc, argv);

  // Q: bcarnes: should run return an arg for the program on exit?
  return 0;
}
