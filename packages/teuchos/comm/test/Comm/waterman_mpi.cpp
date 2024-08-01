// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// MPI-only (no Teuchos::Comm) version of waterman_teuchoscomm.cpp
// This version runs correctly on waterman.
#include <iostream>
#include <fstream>
#include <mpi.h>

void runTest(int ndoubles, MPI_Comm &comm)
{
  int me, np;
  MPI_Comm_rank(comm, &me);
  MPI_Comm_size(comm, &np);

  double *sendBuf = new double[ndoubles];
  for (int i = 0; i < ndoubles; i++) sendBuf[i] = 0.;
  int nMy = ndoubles / np + (me < (ndoubles % np));
  int myBegin = (ndoubles / np) * me;
  myBegin += ((ndoubles % np) < me ? (ndoubles % np) : me);
  int myEnd = myBegin + nMy;
  for (int i = myBegin; i < myEnd; ++i) sendBuf[i] = me;

  double *recvBuf = new double[ndoubles];

  if (me == 0) 
    std::cout << "Trying reduceAll with ndoubles = " << ndoubles << std::endl;

  MPI_Allreduce(sendBuf, recvBuf, ndoubles, MPI_DOUBLE, MPI_SUM, comm);

  delete [] recvBuf;
  delete [] sendBuf;
}


int main(int narg, char **arg)
{
  MPI_Init(&narg, &arg);
  MPI_Comm comm = MPI_COMM_WORLD;

  int me;
  MPI_Comm_rank(comm, &me);

  if (narg != 2) {
    if (me == 0)
      std::cout << "Usage:  a.out [Y|N] \n"
                << "        a.out Y ==> duplicate communicator \n"
                << "        a.out N ==> do not duplicate communicator \n"
                << std::endl;
    return -1;
  }

  bool dupComm = (arg[1][0] == 'Y' ? true : false);
  MPI_Comm commdup;

  if (dupComm) 
    MPI_Comm_dup(comm, &commdup);

  runTest(512, comm);

  if (me == 0) 
    std::cout << "PASSED with " 
              << (dupComm ? "comm duplication " : "no comm duplication ")
              << std::endl;

  if (dupComm) MPI_Comm_free(&commdup);
  MPI_Finalize();
  return 0;
}
