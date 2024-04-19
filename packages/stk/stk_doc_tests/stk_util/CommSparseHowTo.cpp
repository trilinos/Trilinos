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
#include "stk_util/parallel/CommSparse.hpp"    // for CommSparse
#include "stk_util/parallel/Parallel.hpp"      // for MPI_COMM_WORLD, MPI_Comm, ompi_communicator_t
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer
#include "stk_util/stk_config.h"               // for STK_HAS_MPI

#if defined ( STK_HAS_MPI )

//BEGINCommSparse
TEST(ParallelComm, HowToCommunicateOneValue_PackAndCommunicate)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::CommSparse commSparse(comm);

  int myProcId = commSparse.parallel_rank();
  int numProcs = commSparse.parallel_size();

  double sendSomeNumber = 100-myProcId;

  stk::pack_and_communicate(commSparse, [&commSparse, &sendSomeNumber, &myProcId, &numProcs]() {
    for (int proc=0;proc<numProcs;proc++) {
      if ( proc != myProcId ) {
        stk::CommBuffer& proc_buff = commSparse.send_buffer(proc);
        proc_buff.pack<double>(sendSomeNumber);
      }
    }
  });

  for (int proc=0;proc<numProcs;proc++) {
    if ( proc != myProcId ) {
      stk::CommBuffer& dataReceived = commSparse.recv_buffer(proc);
      double val = -1;
      dataReceived.unpack(val);
      EXPECT_EQ(100-proc, val);
    }
  }
}

TEST(ParallelComm, HowToCommunicateOneValue_RawCommSparse)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::CommSparse commSparse(comm);

  int myProcId = commSparse.parallel_rank();
  int numProcs = commSparse.parallel_size();

  double sendSomeNumber = 100-myProcId;

  for(int phase = 0; phase < 2; ++phase) {
    for (int proc=0;proc<numProcs;proc++) {
      if ( proc != myProcId ) {
        stk::CommBuffer& proc_buff = commSparse.send_buffer(proc);
        proc_buff.pack<double>(sendSomeNumber);
      }
    }
    if(phase == 0) {
      commSparse.allocate_buffers();
    }
    else {
      commSparse.communicate();
    }
  }

  for (int proc=0;proc<numProcs;proc++) {
    if ( proc != myProcId ) {
      stk::CommBuffer& dataReceived = commSparse.recv_buffer(proc);
      double val = -1;
      dataReceived.unpack(val);
      EXPECT_EQ(100-proc, val);
    }
  }
}
//ENDCommSparse
//BEGINCommSparse2
TEST(ParallelComm, HowToCommunicateAnArbitraryNumberOfValues)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::CommSparse commSparse(comm);

  int myProcId = commSparse.parallel_rank();
  int numProcs = commSparse.parallel_size();

  double sendSomeNumber = 100-myProcId;

  stk::pack_and_communicate(commSparse, [&commSparse, &sendSomeNumber, &myProcId, &numProcs]() {
    for (int proc=0;proc<numProcs;proc++) {
      if ( proc != myProcId ) {
        stk::CommBuffer& proc_buff = commSparse.send_buffer(proc);
        for (int i=0;i<myProcId;i++) {
          proc_buff.pack<double>(sendSomeNumber+i);
        }
      }
    }
  });

  for (int sourceProc=0; sourceProc < numProcs; sourceProc++) {
    if ( sourceProc != myProcId ) {
      stk::CommBuffer& dataReceived = commSparse.recv_buffer(sourceProc);
      int numItemsReceived = 0;
      while ( dataReceived.remaining() ) {
        double val = -1;
        dataReceived.unpack(val);
        EXPECT_EQ(100-sourceProc+numItemsReceived, val);
        numItemsReceived++;
      }
      int goldNumItemsReceived = sourceProc;
      EXPECT_EQ(goldNumItemsReceived, numItemsReceived);
    }
  }
}
//ENDCommSparse2
#endif
