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
#include <stk_util/parallel/MPITagManager.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <memory>

namespace {

class MPITagManagerForTesting : public stk::MPITagManager
{
  public:
    MPITagManagerForTesting(int deletionGroupSize, int delayCount) :
      MPITagManager(deletionGroupSize, delayCount)
    {}
};

std::shared_ptr<stk::MPITagManager> make_tag_manager_for_testing(int deletionGroupSize, int delayCount)
{
  return std::make_shared<MPITagManagerForTesting>(deletionGroupSize, delayCount);
}

void perf_test(int deletionGroupSize, int delayCount, int ntags, int freeEvery)
{
  auto tagManager = make_tag_manager_for_testing(deletionGroupSize, delayCount);

  std::vector<stk::MPITag> tags;
  tags.reserve(ntags);

  MPI_Barrier(MPI_COMM_WORLD);
  auto t_start = MPI_Wtime();
  for (int i=0; i < ntags; ++i) {
    tags.push_back(tagManager->get_tag(MPI_COMM_WORLD, 1));
    if (i % freeEvery == 0) {
      tags.clear();
    }
  }

  auto elapsed_time = MPI_Wtime() - t_start;

  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if (myRank == 0) {
    std::cout << "elapsed time for " << ntags << ", freed every " << freeEvery << " = " << elapsed_time
              << ", time per tag = " << elapsed_time/ntags << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

TEST(MPITagManager, timing)
{
  const unsigned NUM_RUNS = 5;

  perf_test(32, 16, 256000, 16384);

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  for(unsigned r=0; r<NUM_RUNS; ++r) {
    batchTimer.start_batch_timer();

    perf_test(32, 16, 1024,  1024);
    perf_test(32, 16, 1024,  96);
    perf_test(32, 16, 2048,  2048);
    perf_test(32, 16, 2048,  96);
    perf_test(32, 16, 8192,  8192);
    perf_test(32, 16, 8192,  96);
    perf_test(32, 16, 32768, 16384);
    perf_test(32, 16, 32768, 96);
    perf_test(32, 16, 256000, 16384);
    perf_test(32, 16, 256000, 96);

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_RUNS);
}

} // namespace anonymous

