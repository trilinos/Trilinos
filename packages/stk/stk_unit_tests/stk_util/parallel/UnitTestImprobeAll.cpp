#include "mpi.h"
#include <vector>
#include <cassert>
#include <stdexcept>
#include <iostream>

// comment this out to build this file as its own executable
#define BUILDING_IN_STK

#ifdef BUILDING_IN_STK
#include "gtest/gtest.h"
#endif

namespace {

using DType = int;
const int NVals = 4096;  // use fixed number of elements to simplify this
                         // example.  In the real code the number of
                         // elements is unknown, which is why we use
                         // MPI_Improbe

// compute a unique-ish value (given the receiving rank and the index
// in the buffer, we can always compute who the sender is)
DType getValue(DType rank_src, DType rank_dest, DType idx)
{
  return rank_src + 2*rank_dest + 3*idx;
}

std::vector<std::vector<DType>> getSendData(MPI_Comm comm)
{
  int comm_size, comm_rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  
  std::vector<std::vector<DType>> send_bufs(comm_size);
  for (int dest=0; dest < comm_size; ++dest)
  {
    send_bufs[dest].resize(NVals);
    for (int j=0; j < NVals; ++j)
      send_bufs[dest][j] = getValue(comm_rank, dest, j);
  }

  return send_bufs;
}

void checkRecvData(MPI_Comm comm, std::vector<std::vector<DType>>& recv_bufs)
{
  int comm_size, comm_rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
 
  for (int src=0; src < comm_size; ++src)
    for (int j=0; j < NVals; ++j)
      if (recv_bufs[src][j] != getValue(src, comm_rank, j))
        throw std::runtime_error("error in received data");
}

void doCommunication(MPI_Comm comm)
{
  MPI_Barrier(comm);
  int comm_size, comm_rank, tag=7;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
 

  std::vector<std::vector<DType>> send_bufs = getSendData(comm);
  std::vector<std::vector<DType>> recv_bufs(comm_size);
  std::vector<MPI_Request> send_reqs(comm_size), recv_reqs(comm_size);

  // launch sends
  for (int dest=0; dest < comm_size; ++dest)
    MPI_Isend(send_bufs[dest].data(), NVals*sizeof(DType), MPI_BYTE, dest, tag, comm, &(send_reqs[dest]));


  // using MPI_Irecv works
  /*
  for (int src=0; src < comm_size; ++src)
  {
    recv_bufs[src].resize(NVals);
    MPI_Irecv(recv_bufs[src].data(), NVals*sizeof(DType), MPI_BYTE, src, tag, comm, &(recv_reqs[src]));
  }
  */

  // using MPI_Improbe does not work 

  // Probe for receives
  int nreceived = 0;
  while (true)
  {
    int flag=false;
    MPI_Message message;
    MPI_Status status;
    MPI_Improbe(MPI_ANY_SOURCE, tag, comm, &flag, &message, &status);

    if (flag)
    {
      int src=status.MPI_SOURCE, count=0;
      MPI_Get_count(&status, MPI_BYTE, &count);
      assert(count % sizeof(DType) == 0);
      assert(count == NVals * sizeof(DType));
      //std::cout << "receiving " << count << " bytes from rank " << src << std::endl;

      recv_bufs[src].resize(count / sizeof(DType));
      MPI_Imrecv(recv_bufs[src].data(), count, MPI_BYTE, &message, &(recv_reqs[src]));
      nreceived++;
    }

    if (nreceived == comm_size)
      break;
  }

  MPI_Waitall(send_reqs.size(), send_reqs.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_reqs.size(), recv_reqs.data(), MPI_STATUSES_IGNORE);
  MPI_Barrier(comm);
  checkRecvData(comm, recv_bufs);
}

}  // namespace


#ifdef BUILDING_IN_STK
TEST(Improbe, Improbe_All)
{
  doCommunication(MPI_COMM_WORLD);
}

#else

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  
  doCommunication(MPI_COMM_WORLD);

  MPI_Finalize();
}
#endif
