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

#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_unit_test_utils/timer.hpp>

namespace {

TEST(ParallelVectorConcat, Timing)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  constexpr unsigned numIters = 3;
  constexpr unsigned numRuns = 5;

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  for(unsigned run=0; run<numRuns; ++run) {
    batchTimer.start_batch_timer();

    using Item = std::pair<double,double>;
    const size_t bufSize = size_t(std::numeric_limits<int>::max())/sizeof(Item) + 1000;
    std::vector<Item> localBuf(bufSize);
    std::vector<Item>globalBuf;

    const int myrank = stk::parallel_machine_rank(MPI_COMM_WORLD);
    for (size_t i=0; i < bufSize; ++i) {
      localBuf[i] = std::make_pair(static_cast<double>(i + myrank),static_cast<double>(i + myrank));
    }

    for(unsigned iter=0; iter<numIters; ++iter) {
      stk::parallel_vector_concat(MPI_COMM_WORLD, localBuf, globalBuf);
      EXPECT_EQ(globalBuf.size(), 2*bufSize);

      for (size_t j=0; j < bufSize; ++j) {
        EXPECT_EQ(globalBuf[j], Item(double(j), double(j)));
        EXPECT_EQ(globalBuf[bufSize + j], Item(double(1 + j), double(1 + j)));
      }
    }

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIters);
}

}

