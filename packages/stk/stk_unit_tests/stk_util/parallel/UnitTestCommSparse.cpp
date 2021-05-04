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
#include "stk_util/parallel/CommSparse.hpp"    // for CommSparse, comm_recv_msg_sizes, comm_recv...
#include "stk_util/parallel/Parallel.hpp"      // for parallel_machine_rank, parallel_machine_size
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer
#include "stk_util/stk_config.h"               // for STK_HAS_MPI
#include "stk_util/util/ReportHandler.hpp"     // for ThrowRequireMsg
#include <memory>                              // for allocator_traits<>::value_type
#include <ostream>                             // for basic_ostream::operator<<, operator<<, bas...
#include <vector>                              // for vector

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
      ThrowRequireMsg( recv_bufs[p].size() == expected, "proc "<<myProc<<", recv_bufs["<<p<<"].size()="<<recv_bufs[p].size()<<std::endl);
      EXPECT_EQ(expected, recv_bufs[p].size());
    }
  }
}

TEST(ParallelComm, CommSparse_pair_with_string)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 2) { return; }
  int myRank = stk::parallel_machine_rank(comm);
  int otherRank = 1 - myRank;

  stk::CommSparse commSparse(comm);


  const bool needToUnpackRecvdMessage =
    stk::pack_and_communicate(commSparse, [&commSparse, &myRank, &otherRank]() {
      std::string str = "message from "+std::to_string(myRank)+" to "+std::to_string(otherRank);
      std::pair<std::string,int> pairToSend(str,myRank);
      commSparse.send_buffer(otherRank).pack(pairToSend);
    });

  EXPECT_TRUE(needToUnpackRecvdMessage);

  stk::CommBuffer& buf = commSparse.recv_buffer(otherRank);
  EXPECT_TRUE(buf.remaining() > 0);

  std::pair<std::string,int> expected("message from "+std::to_string(otherRank)+" to "+std::to_string(myRank),otherRank);
  std::pair<std::string,int> recvdPair;
  buf.unpack(recvdPair);
  EXPECT_EQ(expected, recvdPair);
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

TEST(ParallelComm, CommSparse_pack_and_communicate)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  stk::CommSparse commSparse(comm);
  const int numProcs = commSparse.parallel_size();
  if (numProcs == 1) { GTEST_SKIP(); }
  const int myProc = commSparse.parallel_rank();
  const int destProc = myProc+1 == numProcs ? 0 : myProc+1;
  const int srcProc = myProc-1 < 0 ? numProcs-1 : myProc-1;

  stk::pack_and_communicate(commSparse, [&commSparse, &destProc, &myProc]() {
    commSparse.send_buffer(destProc).pack(myProc);
  });

  stk::unpack_communications(commSparse, [&commSparse, &srcProc](int fromProc) {
    if (fromProc == srcProc) {
      int recvData = -1;
      commSparse.recv_buffer(fromProc).unpack(recvData);
      EXPECT_EQ(srcProc, recvData);
    }
  });
}

TEST(ParallelComm, CommSparse_communicate_with_unpack)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  stk::CommSparse commSparse(comm);
  const int numProcs = commSparse.parallel_size();
  if (numProcs == 1) { GTEST_SKIP(); }
  const int myProc = commSparse.parallel_rank();
  const int destProc = myProc+1 == numProcs ? 0 : myProc+1;
  const int srcProc = myProc-1 < 0 ? numProcs-1 : myProc-1;

  commSparse.send_buffer(destProc).pack(myProc);
  commSparse.allocate_buffers();
  commSparse.send_buffer(destProc).pack(myProc);

  commSparse.communicate_with_unpack([&srcProc](int fromProc, stk::CommBuffer& buf) {
    if (fromProc == srcProc) {
      int recvData = -1;
      buf.unpack(recvData);
      EXPECT_EQ(srcProc, recvData);
    }
  });
}

#endif
