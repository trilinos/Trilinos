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
#include <stk_util/util/MemoryTracking.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <cstdlib>

TEST(EntityCommDatabase, TimingInsertErase )
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  std::cout << "HWM at startup: "<<stk::get_max_hwm_across_procs(MPI_COMM_WORLD)<<std::endl;

  const int arbitrarySeed = 42121;
  std::srand(arbitrarySeed);

  const stk::mesh::EntityId minId = 1;
  const stk::mesh::EntityId maxId = 75000;

  auto randomId = [&](){ return minId + std::rand()%(maxId-minId); };

  size_t curMemStart = 0, hwm = 0;
  stk::get_memory_usage(curMemStart, hwm);
  std::cout<<"Start cur-mem: "<<curMemStart<<std::endl;

#ifndef NDEBUG
  const unsigned numIters = 100000;
#else
  const unsigned numIters = 8000000;
#endif
  const unsigned NUM_RUNS = 5;

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    stk::mesh::EntityCommDatabase ecd;
  
    const int owner = 1;
    const int sharer = 2;
    for(unsigned i=0; i<numIters; ++i) {
      stk::mesh::EntityKey key(stk::topology::NODE_RANK, randomId());
      ecd.insert(key, stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, sharer), owner);

      stk::mesh::EntityKey rmKey(stk::topology::NODE_RANK, randomId());
      ecd.erase(rmKey, stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, sharer));
    }

    size_t curMemEnd = 0;
    stk::get_memory_usage(curMemEnd, hwm);
    size_t curMemUsed = curMemEnd - curMemStart;
    std::cout << "EntityCommDatabase Memory Usage: " << stk::human_bytes(curMemUsed) << std::endl;
    std::cout << "EntityCommDatabase (nkeys="<<ecd.num_comm_keys()<<") uses approx "<<(curMemUsed/ecd.num_comm_keys())<<" bytes per key"<<std::endl;

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIters);
}

