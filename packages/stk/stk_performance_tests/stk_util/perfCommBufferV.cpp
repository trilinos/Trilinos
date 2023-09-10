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

#include <string>
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/CommBufferV.hpp>
#include <stk_unit_test_utils/timer.hpp>

TEST(PerfCommBufferV, pack_unpack)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for(unsigned r=0; r<NUM_RUNS; ++r) {
    batchTimer.start_batch_timer();

#ifndef NDEBUG
    const int n = 10000;
#else
    const int n = 1000000;
#endif
    std::vector<stk::CommBufferV> buffers(n);
  
    const int npacks = 100;
    int int1 = 9;
    double double1 = 99.9;
  
    for(int i=0; i<n; ++i) {
      stk::CommBufferV& buf = buffers[i];
      for(int j=0; j<npacks; ++j) {
        buf.pack<double>(double1);
        buf.pack<int>(int1);
        buf.pack<double>(double1);
        buf.pack<int>(int1);
      }
    }
  
    for(int i=0; i<n; ++i) {
      stk::CommBufferV& buf = buffers[i];
      for(int j=0; j<npacks; ++j) {
        buf.unpack<double>(double1);
        buf.unpack<int>(int1);
        buf.unpack<double>(double1);
        buf.unpack<int>(int1);
      }
    }

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_RUNS);
}

