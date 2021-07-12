// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "gtest/gtest.h"
#include "stk_util/parallel/CommBufferV.hpp"    // for CommBufferV
#include "stk_util/parallel/CommNeighbors.hpp"  // for CommNeighbors, STK_MPI_SUPPORTS_NEIGHBOR_...
#include "stk_util/parallel/Parallel.hpp"       // for parallel_machine_rank, parallel_machine_size
#include "stk_util/stk_config.h"                // for STK_HAS_MPI
#include <cstddef>                              // for size_t
#include <vector>                               // for vector

#if defined ( STK_HAS_MPI )

std::vector<int> get_neighbor_procs(int numAllProcs, int localProc, int numNeighbors)
{
    std::vector<int> neighbors;
    for(int proc=localProc-numNeighbors/2; proc<=localProc+numNeighbors/2; ++proc) {
        if (proc >= 0 && proc != localProc && proc < numAllProcs) {
            neighbors.push_back(proc);
        }
    }
    return neighbors;
}

TEST(Parallel, CommNeighborsOneSided_Raw_MPI)
{
#ifdef OMPI_MAJOR_VERSION
#if OMPI_MAJOR_VERSION < 2
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif
#endif

#if !defined(STK_MPI_SUPPORTS_NEIGHBOR_COMM)
  return;
#else

  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = 1;
  MPI_Comm_size(comm, &numProcs);
  if (numProcs != 2) {
    return;
  }
 
  //this test will send from proc 0 to proc 1.
  //proc 0 will not recv anything and proc 1 will not send anything.
  //'sendProcs' means 'procs that I will send to'
  //'recvprocs' means 'procs that I will recv from'
 
  int localProc = 0;
  MPI_Comm_rank(comm, &localProc);
  int otherProc = 1 - localProc;
  std::vector<int> neighbors = {otherProc};
  std::vector<int> sendProcs = {otherProc};
  std::vector<int> recvProcs = {otherProc};
  if (localProc==0) {
    recvProcs.clear();
  }
  else {
    sendProcs.clear();
  }

  MPI_Info info;
  MPI_Info_create(&info);
  int reorder = 0;
  const int* weights = (int*)MPI_UNWEIGHTED;

  MPI_Comm neighborComm;
  MPI_Dist_graph_create_adjacent(comm,
              neighbors.size(), neighbors.data(), weights,
              neighbors.size(), neighbors.data(), weights,
              info, reorder, &neighborComm);
  MPI_Info_free(&info);

  int numItems = 5;
  std::vector<double> sendBuf(numItems, 0), recvBuf;
  std::vector<int> sendCounts(neighbors.size(), 0), recvCounts(neighbors.size());
  std::vector<int> sendDispls(neighbors.size(), 0), recvDispls(neighbors.size());

  for(size_t i=0; i<sendProcs.size(); ++i) {
    sendCounts[i] = numItems;
  }

  if (localProc == 0) {
    for(int i=0; i<numItems; ++i) {
      sendBuf[i] = (localProc + 1 + i);
    }
  }

  MPI_Neighbor_alltoall((void*)sendCounts.data(), 1, MPI_INT,
                        (void*)recvCounts.data(), 1, MPI_INT, neighborComm);

  int totalRecv = 0;
  recvDispls.resize(recvCounts.size());
  for(size_t i=0; i<recvCounts.size(); ++i) {
    recvDispls[i] = totalRecv;
    totalRecv += recvCounts[i];
  }
  recvBuf.resize(totalRecv);

  MPI_Neighbor_alltoallv(
      (void*)sendBuf.data(), sendCounts.data(), sendDispls.data(), MPI_DOUBLE,
      (void*)recvBuf.data(), recvCounts.data(), recvDispls.data(), MPI_DOUBLE, neighborComm);

  for(unsigned i=0; i<recvProcs.size(); ++i) {
      EXPECT_EQ(1, localProc);//should only recv on proc 1

      const double* recvdFromProcI = recvBuf.data() + recvDispls[i];
      EXPECT_EQ(numItems, recvCounts[i]);

      for(int item=0; item<numItems; ++item) {
          double val = recvdFromProcI[item];
          double expected = recvProcs[i] + 1 + item;
          EXPECT_EQ(expected, val);
      }
  }
#endif
}

TEST(Parallel, CommNeighborsOneSided)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs != 2) {
    return;
  }
 
  int localProc = stk::parallel_machine_rank(comm);
  int otherProc = 1 - localProc;
  std::vector<int> sendProcs = {otherProc};
  std::vector<int> recvProcs = {otherProc};
  if (localProc==0) {
    recvProcs.clear();
  }
  else {
    sendProcs.clear();
  }

  stk::CommNeighbors commNeighbors(comm, sendProcs, recvProcs);

  int numItems = 5;
  for(int proc : sendProcs) {
    stk::CommBufferV& buf = commNeighbors.send_buffer(proc);
    for(int i=0; i<numItems; ++i) {
      buf.pack<double>(localProc+1+i);
    }
  }

  commNeighbors.communicate();

  for(int proc : recvProcs) {
      stk::CommBufferV& buf = commNeighbors.recv_buffer(proc);

      int remaining = buf.size_in_bytes()/sizeof(double);
      EXPECT_EQ(numItems, remaining);

      for(int i=0; i<numItems; ++i) {
          double val = 0;
          buf.unpack<double>(val);
          double expected = proc+1+i;
          EXPECT_EQ(expected, val);
      }
  }
}

TEST(Parallel, CommNeighbors)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs == 1) {
    return;
  }
 
  int myProc = stk::parallel_machine_rank(comm);
  int numNeighbors = 5;
  std::vector<int> neighborProcs = get_neighbor_procs(numProcs, myProc, numNeighbors);
  int numItems = 20;
 
  stk::parallel_machine_barrier(comm);
 
  stk::CommNeighbors commNeighbors(comm, neighborProcs);

  EXPECT_EQ(numProcs, commNeighbors.parallel_size());
  EXPECT_EQ(myProc, commNeighbors.parallel_rank());

  for(int proc : neighborProcs) {
      stk::CommBufferV& proc_buff = commNeighbors.send_buffer(proc);
      for(int i=0; i<numItems; ++i) {
          proc_buff.pack<double>(myProc+1+i);
      }
  }

  commNeighbors.communicate();

  for(int proc : neighborProcs) {
      stk::CommBufferV& proc_buff = commNeighbors.recv_buffer(proc);

      int remaining = proc_buff.size_in_bytes()/sizeof(double);
      EXPECT_EQ(numItems, remaining);

      for(int i=0; i<numItems; ++i) {
          double val = 0;
          proc_buff.unpack<double>(val);
          double expected = proc+1+i;
          EXPECT_EQ(expected, val);
      }
  }
}

TEST(Parallel, CommNeighborsResetBuffers) {
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  std::vector<int> send, recv;
  int numProcs = stk::parallel_machine_size(comm);
  int rank= stk::parallel_machine_rank(comm);
  if (numProcs == 1) {
    return;
  }
  if(rank == 0) {
    for(int p = 1; p < numProcs; ++p) {
      send.push_back(p);
    }
  } else {
    recv.push_back(0);
  }

  stk::CommNeighbors commNeighbors(comm, send, recv);

  for(int i = 0; i < 3; ++i) {
    commNeighbors.reset_buffers();
    for(int p : commNeighbors.send_procs()) {
      commNeighbors.send_buffer(p).pack<int>(i);
    }

    commNeighbors.communicate();

    for(int p : commNeighbors.recv_procs()) {
      int val = 0;
      commNeighbors.recv_buffer(p).unpack<int>(val);
      EXPECT_EQ(i, val);
      EXPECT_EQ(0u, commNeighbors.recv_buffer(p).size_in_bytes());
    }
  }
}

#endif
