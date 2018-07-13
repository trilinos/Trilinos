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
 
  stk::CommNeighbors commneighbors(comm, neighborProcs);

  for(int proc : neighborProcs) {
      stk::CommBufferV& proc_buff = commneighbors.send_buffer(proc);
      for(int i=0; i<numItems; ++i) {
          proc_buff.pack<double>(myProc+1+i);
      }
  }

  commneighbors.communicate();

  for(int proc : neighborProcs) {
      stk::CommBufferV& proc_buff = commneighbors.recv_buffer(proc);

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
