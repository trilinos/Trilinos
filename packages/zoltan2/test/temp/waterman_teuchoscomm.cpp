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
  memset(sendBuf, 0, sizeof(double) * ndoubles);
  int nMy = ndoubles / np + (me < (ndoubles % np));
  int myBegin = (ndoubles / np) * me;
  myBegin += ((ndoubles % np) < me ? (ndoubles % np) : me);
  int myEnd = myBegin + nMy;
  for (int i = myBegin; i < myEnd; ++i) sendBuf[i] = me;

  double *recvBuf = new double[ndoubles];

  if (me == 0) 
    std::cout << "Trying reduceAll with ndoubles = " << ndoubles << std::endl;

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

  if (narg != 2 && me == 0) {
    std::cout << "Usage:  a.out [0|1] \n"
              << "        a.out 1 ==> duplicate communicator \n"
              << "        a.out 0 ==> do not duplicate communicator \n"
              << std::endl;
    return -1;
  }

  int dupComm = atoi(arg[1]);

  if (dupComm) 
    Teuchos::RCP<const Teuchos::Comm<int> > commdup = comm->duplicate();

  runTest(512, comm);

  if (me == 0) 
    std::cout << "PASSED with " 
              << (dupComm ? "comm duplication " : "no comm duplication ")
              << std::endl;
  return 0;
}
