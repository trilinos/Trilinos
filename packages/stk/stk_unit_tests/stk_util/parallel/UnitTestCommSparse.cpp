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

#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll
#include <stk_util/parallel/CommSparse.hpp>  // for comm_recv_sizes
#include <stk_util/parallel/MPI.hpp>
#include <gtest/gtest.h>
#include <vector>                       // for vector
#include <stk_util/stk_config.h>
#include <limits>

#if defined ( STK_HAS_MPI )

TEST(ParallelComm, comm_recv_msg_sizes)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs == 1) {
    return;
  }
  int myProc = stk::parallel_machine_rank(comm);

  std::vector<stk::CommBuffer> send_bufs(numProcs), recv_bufs(numProcs);
  std::vector<int> send_procs, recv_procs;
  
  for(int p=myProc-2; p<=myProc+2; ++p) {
    if (p >= 0 && p < numProcs) {
      send_bufs[p].set_size(myProc+1);
      send_procs.push_back(p);
      recv_procs.push_back(p);
    }
  }

  stk::comm_recv_msg_sizes(comm, send_procs, recv_procs, send_bufs, recv_bufs);

  for(int p=0; p<numProcs; ++p) {
    if (p < (myProc-2) || p > (myProc+2)) {
      EXPECT_EQ(0u, recv_bufs[p].size());
    }
    else {
      EXPECT_EQ((unsigned)p+1, recv_bufs[p].size());
    }
  }
}

TEST(ParallelComm, comm_recv_procs_and_msg_sizes)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs == 1) {
    return;
  }
  int myProc = stk::parallel_machine_rank(comm);

  std::vector<stk::CommBuffer> send_bufs(numProcs), recv_bufs(numProcs);
  
  for(int p=myProc-2; p<=myProc+2; ++p) {
    if (p >= 0 && p < numProcs) {
      send_bufs[p].set_size(myProc+1);
    }
  }

  std::vector<int> send_procs, recv_procs;
  stk::comm_recv_procs_and_msg_sizes(comm, send_bufs, recv_bufs, send_procs, recv_procs);

  for(int p=0; p<numProcs; ++p) {
    if (p < (myProc-2) || p > (myProc+2)) {
      EXPECT_EQ(0u, recv_bufs[p].size());
    }
    else {
      unsigned expected = p+1;
      if (recv_bufs[p].size() != expected) {
std::ostringstream msg;
msg<<"proc "<<myProc<<", recv_bufs["<<p<<"].size()="<<recv_bufs[p].size()<<std::endl;
std::cout<<msg.str()<<std::endl;
      }
      EXPECT_EQ(expected, recv_bufs[p].size());
    }
  }
}

TEST(ParallelComm, CommSparse)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs == 1) {
    return;
  }
 
  int myProc = stk::parallel_machine_rank(comm);
  int numOtherProcs = 2;
  int numItems = 20;
 
  stk::parallel_machine_barrier(comm);
 
  stk::CommSparse commsparse(comm);

  for(int phase = 0; phase < 2; ++phase) {
    for(int proc=myProc-numOtherProcs/2; proc<=myProc+numOtherProcs/2; ++proc) {
      if (proc >= 0 && proc != myProc && proc < numProcs) {
        stk::CommBuffer& proc_buff = commsparse.send_buffer(proc);
        for(int i=0; i<numItems; ++i) {
          proc_buff.pack<double>(myProc+1+i);
        }
      }
    }

    if (phase == 0) {
      commsparse.allocate_buffers();
    }
    else {
      commsparse.communicate();
    }
  }

  for(int proc=0; proc<numProcs; ++proc) {
    if (proc != myProc && (proc >= myProc-numOtherProcs/2) && (proc <= myProc+numOtherProcs/2)) {
      stk::CommBuffer& proc_buff = commsparse.recv_buffer(proc);

      int remaining = proc_buff.remaining()/sizeof(double);
      EXPECT_EQ(numItems, remaining);

      for(int i=0; i<numItems; ++i) {
        double val = 0;
        proc_buff.unpack<double>(val);
        double expected = proc+1+i;
        EXPECT_EQ(expected, val);
      }
    }
    else {
      stk::CommBuffer& proc_buff = commsparse.recv_buffer(proc);

      int remaining = proc_buff.remaining();
      EXPECT_EQ(0, remaining);
    }
  }
}

TEST(ParallelComm, CommSparse_set_procs)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs == 1) {
    return;
  }
 
  int myProc = stk::parallel_machine_rank(comm);
  int numOtherProcs = 2;
  int numItems = 20;
 
  stk::parallel_machine_barrier(comm);
 
  stk::CommSparse commsparse(comm);

  std::vector<int> send_procs, recv_procs;

  for(int phase = 0; phase < 2; ++phase) {
    for(int proc=myProc-numOtherProcs/2; proc<=myProc+numOtherProcs/2; ++proc) {
      if (proc >= 0 && proc != myProc && proc < numProcs) {
        stk::CommBuffer& proc_buff = commsparse.send_buffer(proc);
        send_procs.push_back(proc);
        recv_procs.push_back(proc);
        for(int i=0; i<numItems; ++i) {
          proc_buff.pack<double>(myProc+1+i);
        }
      }
    }

    if (phase == 0) {
      commsparse.allocate_buffers(send_procs, recv_procs);
    }
    else {
      commsparse.communicate();
    }
  }

  for(int proc=0; proc<numProcs; ++proc) {
    if (proc != myProc && (proc >= myProc-numOtherProcs/2) && (proc <= myProc+numOtherProcs/2)) {
      stk::CommBuffer& proc_buff = commsparse.recv_buffer(proc);

      int remaining = proc_buff.remaining()/sizeof(double);
      EXPECT_EQ(numItems, remaining);

      for(int i=0; i<numItems; ++i) {
        double val = 0;
        proc_buff.unpack<double>(val);
        double expected = proc+1+i;
        EXPECT_EQ(expected, val);
      }
    }
    else {
      stk::CommBuffer& proc_buff = commsparse.recv_buffer(proc);

      int remaining = proc_buff.remaining();
      EXPECT_EQ(0, remaining);
    }
  }
}

#endif
