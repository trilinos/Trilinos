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

// send only to self
void doCommunication(MPI_Comm comm)
{
  MPI_Barrier(comm);
  int comm_size, comm_rank, tag=7;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
 
  std::vector<DType> send_buf(NVals), recv_buf(NVals);
  MPI_Request send_req, recv_req;

  for (int i=0; i < NVals; ++i)
    send_buf[i] = i + comm_rank;


  // post send
  MPI_Isend(send_buf.data(), NVals*sizeof(DType), MPI_BYTE, comm_rank, tag, comm, &send_req);

  // Probe for receives
  while (true)
  {
    int flag=false;
    MPI_Message message;
    MPI_Status status;
    MPI_Improbe(MPI_ANY_SOURCE, tag, comm, &flag, &message, &status);

    if (flag)
    {
      int count=0;
      MPI_Get_count(&status, MPI_BYTE, &count);
      assert(count == NVals * sizeof(DType));
      assert(status.MPI_SOURCE == comm_rank);

      MPI_Imrecv(recv_buf.data(), count, MPI_BYTE, &message, &recv_req);
      break;
    }
  }

  MPI_Wait(&send_req, MPI_STATUS_IGNORE);
  MPI_Wait(&recv_req, MPI_STATUS_IGNORE);

  // check result
  for (int i=0; i < NVals; ++i)
    if (send_buf[i] != recv_buf[i])
      throw std::runtime_error("recv buffer value not correct");

  MPI_Barrier(comm);
}

}

#ifdef BUILDING_IN_STK
TEST(Improbe, Improbe_self)
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
