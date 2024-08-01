// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test to exercise fix to Teuchos::Comm's createSubcommunicator, which
// had leaked memory.
// The fix added MPI_Comm_free to the opaqueWrappers. 
// 8/2018  This test hangs of platform waterman.

#include <stdio.h>
#include <mpi.h>
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"


int main(int narg, char **arg)
{
  Teuchos::GlobalMPISession mpiSession(&narg,&arg);

  Teuchos::RCP<const Teuchos::Comm<int> >
    comm = Teuchos::DefaultComm<int>::getComm();
  int me = comm->getRank();
  int np = comm->getSize();

  int niter = 4;
  int *ids = new int[np/2+1];
  for (int i = 0; i < np/2+1; i++) ids[i] = i;
  Teuchos::ArrayView<const int> list(ids, np/2+1);

  for (int i = 0; i < niter; i++) {
    Teuchos::RCP<const Teuchos::Comm<int> > a 
                                            = comm->createSubcommunicator(list);
    printf("iteration %d -- weak: %d  strong: %d total: %d\n",
            i, a.weak_count(), a.strong_count(), a.total_count());
  }
  delete [] ids;
  if (me == 0)
    printf("\nPASS\n");

  return 0;
}
