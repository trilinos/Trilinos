#include <iostream>
#include <fstream>
#include <mpi.h>

void runTest(int ndoubles, MPI_Comm &comm)
{
  int me, np;
  MPI_Comm_rank(comm, &me);
  MPI_Comm_size(comm, &np);

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

  if (narg != 2 && me == 0) {
    std::cout << "Usage:  a.out [0|1] \n"
              << "        a.out 1 ==> duplicate communicator \n"
              << "        a.out 0 ==> do not duplicate communicator \n"
              << std::endl;
    return -1;
  }

  int dupComm = atoi(arg[1]);
  MPI_Comm commdup;

  if (dupComm) 
    MPI_Comm_dup(comm, &commdup);

  runTest(512, comm);

  if (me == 0) 
    std::cout << "PASSED with " 
              << (dupComm ? "comm duplication " : "no comm duplication ")
              << std::endl;

  MPI_Finalize();
  return 0;
}
