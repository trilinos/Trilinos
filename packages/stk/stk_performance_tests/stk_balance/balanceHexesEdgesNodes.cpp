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
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_balance/balance.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/timer.hpp>

class BalanceHexesEdgesNodes : public stk::unit_test_util::MeshFixture
{
public:
  BalanceHexesEdgesNodes()
    : stk::unit_test_util::MeshFixture()
  { }

  void setup_host_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
#ifdef NDEBUG
    setup_mesh("generated:200x200x64", auraOption);
#else
    setup_mesh("generated:10x10x100", auraOption);
#endif

    stk::mesh::create_edges(get_bulk());
  }

  void perform_balance(const std::string& decompMethod)
  {
    stk::balance::GraphCreationSettings balanceSettings;
    balanceSettings.setDecompMethod(decompMethod);

    stk::balance::balanceStkMesh(balanceSettings, get_bulk());

    balanceSettings.setUseNodeBalancer(true);
    balanceSettings.setNodeBalancerTargetLoadBalance(1);
    balanceSettings.setNodeBalancerMaxIterations(4);

    stk::balance::balanceStkMeshNodes(balanceSettings, get_bulk());
  }

  void run_balance_performance_test(stk::unit_test_util::BatchTimer& batchTimer,
                                    unsigned NUM_RUNS,
                                    const std::string& decompMethod,
                                    stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    for (unsigned j = 0; j < NUM_RUNS; j++) {
      setup_host_mesh(auraOption);

      batchTimer.start_batch_timer();
    
      perform_balance(decompMethod);

      batchTimer.stop_batch_timer();

      reset_mesh();
    }
  }
};

TEST_F(BalanceHexesEdgesNodes, RIB_withAura)
{
  if (get_parallel_size() != 8) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();

  run_balance_performance_test(batchTimer, NUM_RUNS, "rib", stk::mesh::BulkData::AUTO_AURA);

  batchTimer.print_batch_timing(NUM_RUNS);
}

TEST_F(BalanceHexesEdgesNodes, RIB_noAura)
{
  if (get_parallel_size() != 8) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();

  run_balance_performance_test(batchTimer, NUM_RUNS, "rib", stk::mesh::BulkData::NO_AUTO_AURA);

  batchTimer.print_batch_timing(NUM_RUNS);
}

TEST_F(BalanceHexesEdgesNodes, parmetis_withAura)
{
  if (get_parallel_size() != 8) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();

  run_balance_performance_test(batchTimer, NUM_RUNS, "parmetis", stk::mesh::BulkData::AUTO_AURA);

  batchTimer.print_batch_timing(NUM_RUNS);
}

TEST_F(BalanceHexesEdgesNodes, parmetis_noAura)
{
  if (get_parallel_size() != 8) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();

  run_balance_performance_test(batchTimer, NUM_RUNS, "parmetis", stk::mesh::BulkData::NO_AUTO_AURA);

  batchTimer.print_batch_timing(NUM_RUNS);
}

