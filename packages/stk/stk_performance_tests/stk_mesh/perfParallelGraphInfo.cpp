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
#include <stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp>
#include <stk_performance_tests/stk_mesh/timer.hpp>
#include <cstdlib>

TEST(ParallelGraphInfo, TimingInsertErase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  std::cout << "HWM at startup: "<<stk::get_max_hwm_across_procs(MPI_COMM_WORLD)<<std::endl;
  std::cout << "sizeof(GraphEdge): "<<sizeof(stk::mesh::GraphEdge)<<", sizeof(ParallelInfo): "<<sizeof(stk::mesh::impl::ParallelInfo)<<std::endl;

  const int arbitrarySeed = 42121;
  std::srand(arbitrarySeed);

  const stk::mesh::impl::LocalId minId = 0;
  const stk::mesh::impl::LocalId maxId = 75000;

  auto randomId = [&](){ return minId + std::rand()%(maxId-minId); };

  size_t curMemStart = 0, hwm = 0;
  stk::get_memory_usage(curMemStart, hwm);
  std::cout<<"Start cur-mem: "<<curMemStart<<std::endl;

#ifndef NDEBUG
  const unsigned numIters = 100000;
#else
  const unsigned numIters = 5000000;
#endif
  const unsigned NUM_RUNS = 5;

  stk::performance_tests::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    const int proc = 1;
    stk::mesh::ParallelInfoForGraphEdges pllGraphInfo(proc);
  
    const int otherProc = 2;

    for(unsigned i=0; i<numIters; ++i) {
      {
        stk::mesh::impl::LocalId localElemId = randomId();
        int localSide = 1;
        stk::mesh::impl::LocalId negativeRemoteId = -localElemId;
        int remoteSide = 4;
        stk::mesh::GraphEdge graphEdge(localElemId, localSide, negativeRemoteId, remoteSide);

        stk::mesh::impl::ParallelInfo parInfo(otherProc, stk::mesh::INVALID_PERMUTATION,
                                              stk::topology::HEX_8);

        pllGraphInfo.insert_parallel_info_for_graph_edge(graphEdge, parInfo);
      }

      {
        stk::mesh::impl::LocalId localElemId = randomId();
        int localSide = 1;
        stk::mesh::impl::LocalId negativeRemoteId = -localElemId;
        int remoteSide = 4;
        stk::mesh::GraphEdge graphEdge(localElemId, localSide, negativeRemoteId, remoteSide);

        pllGraphInfo.erase_parallel_info_for_graph_edge(graphEdge);
      }
    }

    size_t curMemEnd = 0;
    stk::get_memory_usage(curMemEnd, hwm);
    size_t curMemUsed = curMemEnd - curMemStart;
    std::cout << "ParallelInfoForGraphEdges Memory Usage: " << stk::human_bytes(curMemUsed) << std::endl;
    size_t nedges = pllGraphInfo.get_parallel_graph_info().size();
    std::cout << "ParallelInfoForGraphEdges (nedges="<<nedges<<") uses approx "<<(curMemUsed/nedges)<<" bytes per edge"<<std::endl;

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIters);
}

