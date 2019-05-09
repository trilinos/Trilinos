// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <gtest/gtest.h>
#include <vector>                       // for vector
#include <stk_util/stk_config.h>
#include <limits>
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_util/parallel/CommNeighbors.hpp>
#include <stk_util/parallel/MPI.hpp>

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
#if (defined(OMPI_MAJOR_VERSION) && (OMPI_MAJOR_VERSION < 2 ))
//this test doesn't pass with open-mpi 1.10 but does pass
//with intel-mpi and (presumably) newer open-mpi versions.
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

  //cppreference.com says vector::data() may be nullptr if vector is empty but not guaranteed
  const int* sendProcsPtr = sendProcs.size() > 0 ? sendProcs.data() : nullptr;
  const int* recvProcsPtr = recvProcs.size() > 0 ? recvProcs.data() : nullptr;

  MPI_Comm neighborComm;
  MPI_Dist_graph_create_adjacent(comm,
              recvProcs.size(), recvProcsPtr, weights,
              sendProcs.size(), sendProcsPtr, weights,
              info, reorder, &neighborComm);
  MPI_Info_free(&info);

  int numItems = 5;
  std::vector<double> sendBuf(sendProcs.size()*numItems, 0), recvBuf;
  std::vector<int> sendCounts(sendProcs.size(), numItems), recvCounts(recvProcs.size());
  std::vector<int> sendDispls(sendProcs.size(), 0), recvDispls(recvProcs.size());

  if (localProc == 0) {
    for(int i=0; i<numItems; ++i) {
      sendBuf[i] = (localProc + 1 + i);
    }
  }

  const int* sendCountsPtr = sendCounts.size() > 0 ? sendCounts.data() : nullptr;
  const int* recvCountsPtr = recvCounts.size() > 0 ? recvCounts.data() : nullptr;

  MPI_Neighbor_alltoall((void*)sendCountsPtr, 1, MPI_INT,
                        (void*)recvCountsPtr, 1, MPI_INT, neighborComm);

  int totalRecv = 0;
  recvDispls.resize(recvCounts.size());
  for(size_t i=0; i<recvCounts.size(); ++i) {
    recvDispls[i] = totalRecv;
    totalRecv += recvCounts[i];
  }
  recvBuf.resize(totalRecv);

  const double* sendBufPtr = sendBuf.size() > 0 ? sendBuf.data() : nullptr;
  const double* recvBufPtr = recvBuf.size() > 0 ? recvBuf.data() : nullptr;
  const int* sendDisplsPtr = sendDispls.size() > 0 ? sendDispls.data() : nullptr;
  const int* recvDisplsPtr = recvDispls.size() > 0 ? recvDispls.data() : nullptr;

  MPI_Neighbor_alltoallv(
      (void*)sendBufPtr, sendCountsPtr, sendDisplsPtr, MPI_DOUBLE,
      (void*)recvBufPtr, recvCountsPtr, recvDisplsPtr, MPI_DOUBLE, neighborComm);

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

#endif
