// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Program demonstrating the Teuchos::Comm issues causing #3331.
// On waterman, this test passes if it skips the call to 
// Teuchos::Comm::duplicate, but segfaults if Teuchos::Comm::duplicate 
// is called.
#include <iostream>
#include <fstream>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

void runTest(int ndoubles, Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  int me = comm->getRank();
  int np = comm->getSize();

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

  std::cout << *comm << std::endl;

  Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM,
                                   ndoubles, sendBuf, recvBuf);

  delete [] recvBuf;
  delete [] sendBuf;
}


int main(int narg, char **arg)
{
  Teuchos::GlobalMPISession mpiSession(&narg,&arg);
  Teuchos::RCP<const Teuchos::Comm<int> >
    comm = Teuchos::DefaultComm<int>::getComm();

  int me = comm->getRank();

  if (narg != 2) {
    if (me == 0)
      std::cout << "Usage:  a.out [Y|N] \n"
                << "        a.out Y ==> duplicate communicator \n"
                << "        a.out N ==> do not duplicate communicator \n"
                << std::endl;
    return -1;
  }

  bool dupComm = (arg[1][0] == 'Y' ? true : false);

  if (dupComm) 
    Teuchos::RCP<const Teuchos::Comm<int> > commdup = comm->duplicate();

  runTest(512, comm);

  if (me == 0) 
    std::cout << "PASSED with " 
              << (dupComm ? "comm duplication " : "no comm duplication ")
              << std::endl;
  return 0;
}
