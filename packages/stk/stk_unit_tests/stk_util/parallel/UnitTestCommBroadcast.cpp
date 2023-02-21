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
#include "stk_util/stk_config.h"
#include "stk_util/parallel/Parallel.hpp"      // for parallel_machine_rank, parallel_machine_size
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBroadcast
#include <ostream>                             // for basic_ostream::operator<<, operator<<, bas...
#include <vector>                              // for vector

TEST(ParallelComm, CommBroadcast_string)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 3) { GTEST_SKIP(); }

  const int rootRank = 0;
  stk::CommBroadcast commBroadcast(comm, rootRank);

  EXPECT_EQ(rootRank, commBroadcast.root_rank());

  stk::pack_and_communicate(commBroadcast, [&]() {
    if (commBroadcast.parallel_rank() == rootRank) {
      std::string str = "message from "+std::to_string(commBroadcast.parallel_rank());
      commBroadcast.send_buffer().pack(str);
    }
  });

  stk::CommBuffer& buf = commBroadcast.recv_buffer();
  EXPECT_TRUE(buf.remaining() > 0);

  std::string expected("message from 0");
  std::string recvdStr;
  buf.unpack(recvdStr);
  EXPECT_EQ(expected, recvdStr);
}

