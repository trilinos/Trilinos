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
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_util/parallel/OutputStreams.hpp>
#include <map>
#include <vector>
#include <algorithm>

namespace
{

class BalanceSettingsTester : public stk::balance::GraphCreationSettings
{
public:
  BalanceSettingsTester(const std::string& decompMethod)
    : localMethod(decompMethod) { }
  virtual ~BalanceSettingsTester() = default;

  virtual std::string getDecompMethod() const override { return localMethod; }

private:
  const std::string& localMethod;
};

using CohesiveElementBlocks = std::vector<std::string>;

class TestCohesiveElements : public stk::unit_test_util::MeshFixture
{
protected:

  void set_up_1x1x9_mesh_one_block()
  {
    setup_mesh("generated:1x1x9", stk::mesh::BulkData::AUTO_AURA);
    stk::mesh::Part * block1 = &get_meta().declare_part("block_1");
    move_elements_into_part(block1, {1,2,3,4,5,6,7,8,9});
  }

  void set_up_1x1x9_mesh_two_blocks()
  {
    setup_mesh("generated:1x1x9", stk::mesh::BulkData::AUTO_AURA);
    stk::mesh::Part * block2 = &get_meta().declare_part("block_2");
    move_elements_into_part(block2, {5});
    stk::mesh::Part * block1 = &get_meta().declare_part("block_1");
    move_elements_into_part(block1, {1,2,3,4,6,7,8,9});
  }

  void move_elements_into_part(stk::mesh::Part * addPart, const std::vector<unsigned> & elemIds)
  {
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    stk::mesh::EntityVector elemsToMove;
    for (const unsigned & elemId : elemIds) {
      elemsToMove.push_back(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    get_bulk().batch_change_entity_parts(elemsToMove, {addPart}, {block1});
  }

  void run_decomposition_two_procs(const std::string & decompType,
                                   const CohesiveElementBlocks & cohesiveElements,
                                   const std::vector<unsigned> & expectedElemsOnEachProc)
  {
    BalanceSettingsTester balanceSettings(decompType);
    for (const auto & block : cohesiveElements) {
      balanceSettings.setCohesiveElements(block);
    }

    std::vector<stk::mesh::Selector> selectors = {get_bulk().mesh_meta_data().universal_part()};
    const int numProcs = 2;

    stk::mesh::EntityProcVec decomp;
    stk::set_outputP0(&stk::outputNull());
    stk::balance::internal::calculateGeometricOrGraphBasedDecomp(get_bulk(), selectors,
                                                                 get_bulk().parallel(), numProcs,
                                                                 balanceSettings, decomp);
    stk::reset_default_output_streams();

    check_decomposition(decomp, expectedElemsOnEachProc);
  }

  void check_decomposition(const stk::mesh::EntityProcVec & decomp, const std::vector<unsigned> & expectedElemsOnEachProc)
  {
    std::vector<unsigned> elemCounts(expectedElemsOnEachProc.size(), 0);
    for (const stk::mesh::EntityProc & entityProc : decomp) {
      const int proc = entityProc.second;
      elemCounts[proc]++;
    }
    std::sort(elemCounts.begin(), elemCounts.end());

    std::vector<unsigned> expectedElemCounts = expectedElemsOnEachProc;
    std::sort(expectedElemCounts.begin(), expectedElemCounts.end());

    ASSERT_EQ(elemCounts.size(), expectedElemCounts.size());
    EXPECT_EQ(elemCounts, expectedElemCounts);
  }
};

TEST_F(TestCohesiveElements, OneBlockNoCohesiveElementsParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x9_mesh_one_block();
  const CohesiveElementBlocks cohesiveElements {};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 5};
  run_decomposition_two_procs("parmetis", cohesiveElements, expectedElemsOnEachProc);
}

TEST_F(TestCohesiveElements, OneBlockAllCohesiveElementsParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x9_mesh_one_block();
  const CohesiveElementBlocks cohesiveElements {"block_1"};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 5};
  run_decomposition_two_procs("parmetis", cohesiveElements, expectedElemsOnEachProc);
}

TEST_F(TestCohesiveElements, TwoBlocksNoCohesiveElementsParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x9_mesh_two_blocks();
  const CohesiveElementBlocks cohesiveElements {};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 5};
  run_decomposition_two_procs("parmetis", cohesiveElements, expectedElemsOnEachProc);
}

TEST_F(TestCohesiveElements, TwoBlocksAllCohesiveElementsParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x9_mesh_two_blocks();
  const CohesiveElementBlocks cohesiveElements {"block_1", "block_2"};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 5};
  run_decomposition_two_procs("parmetis", cohesiveElements, expectedElemsOnEachProc);
}

TEST_F(TestCohesiveElements, TwoBlocksCohesiveElementsLayerParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x9_mesh_two_blocks();
  const CohesiveElementBlocks cohesiveElements {"block_2"};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 5};
  run_decomposition_two_procs("parmetis", cohesiveElements, expectedElemsOnEachProc);
}

}
