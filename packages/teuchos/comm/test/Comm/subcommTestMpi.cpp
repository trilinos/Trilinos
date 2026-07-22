// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// MPI-only version of subcommTestTeuchosComm.cpp

#include <iostream>
#include <mpi.h>


int main(int narg, char **arg)
{
  MPI_Init(&narg, &arg);
  MPI_Comm comm = MPI_COMM_WORLD;

  int me, np;
  MPI_Comm_rank(comm, &me);
  MPI_Comm_size(comm, &np);

  int niter = 4;
  int *ids = new int[np/2+1];
  for (int i = 0; i < np/2+1; i++) ids[i] = i;
  
  for (int i = 0; i < niter; i++) {
    MPI_Comm a;
    MPI_Group cgrp, agrp;
    MPI_Comm_group(comm, &cgrp);
    MPI_Group_incl(cgrp, np/2+1, ids, &agrp);

    MPI_Comm_create(comm, agrp, &a);

    MPI_Group_free(&agrp);
    MPI_Group_free(&cgrp);

    if (a != MPI_COMM_NULL) {
      int anp;
      MPI_Comm_size(a, &anp);
      std::cout << me << " Iteration " << i << " New comm has " << anp << " ranks"
                << std::endl;
      MPI_Comm_free(&a);
    }
    else {
      std::cout << me << " not in new communicator" << std::endl;
    }
  }
  delete [] ids;
  if (me == 0)
    std::cout << "PASS" << std::endl;

  MPI_Finalize();
  return 0;
}
