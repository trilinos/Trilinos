#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
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

using BlockWeightsMap = std::map<std::string, double>;

class TestBlockWeights : public stk::unit_test_util::MeshFixture
{
protected:
  void set_up_1x1x8_mesh_one_block()
  {
    setup_mesh("generated:1x1x8", stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void set_up_1x1x8_mesh_two_blocks()
  {
    setup_mesh("generated:1x1x8", stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Part * block2 = &get_meta().declare_part("block_2");
    move_elements_into_part(block2, {5, 6, 7, 8});
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
                                   const BlockWeightsMap & weightForBlock,
                                   const std::vector<unsigned> & expectedElemsOnEachProc)
  {
    BalanceSettingsTester balanceSettings(decompType);
    for (const auto & weights : weightForBlock) {
      balanceSettings.setVertexWeightBlockMultiplier(weights.first, weights.second);
    }

    std::vector<stk::mesh::Selector> selectors = {get_bulk().mesh_meta_data().universal_part()};
    const int numProcs = 2;

    stk::mesh::EntityProcVec decomp;
    testing::internal::CaptureStdout();
    stk::balance::internal::calculateGeometricOrGraphBasedDecomp(get_bulk(), selectors,
                                                                 get_bulk().parallel(), numProcs,
                                                                 balanceSettings, decomp);
    testing::internal::GetCapturedStdout();

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


TEST_F(TestBlockWeights, DefaultWeightsOneBlockParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_one_block();
  const BlockWeightsMap weightForBlock {};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("parmetis", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, DefaultWeightsOneBlockRCB)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_one_block();
  const BlockWeightsMap weightForBlock {};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("rcb", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, DefaultWeightsTwoBlocksParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_two_blocks();
  const BlockWeightsMap weightForBlock {};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("parmetis", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, DefaultWeightsTwoBlocksRCB)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_two_blocks();
  const BlockWeightsMap weightForBlock {};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("rcb", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, MixedDefaultAndSpecifiedWeightsTwoBlocksParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_two_blocks();
  const BlockWeightsMap weightForBlock {{"block_2", 1.0}};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("parmetis", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, MixedDefaultAndSpecifiedWeightsTwoBlocksRCB)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_two_blocks();
  const BlockWeightsMap weightForBlock {{"block_2", 1.0}};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("rcb", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, SpecifiedWeightOneBlockParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_one_block();
  const BlockWeightsMap weightForBlock {{"block_1", 5.0}};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("parmetis", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, SpecifiedWeightOneBlockRCB)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_one_block();
  const BlockWeightsMap weightForBlock {{"block_1", 5.0}};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("rcb", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, UniformWeightsTwoBlocksParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_two_blocks();
  const BlockWeightsMap weightForBlock {{"block_1", 5.0}, {"block_2", 5.0}};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("parmetis", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, UniformWeightsTwoBlocksRCB)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_two_blocks();
  const BlockWeightsMap weightForBlock {{"block_1", 5.0}, {"block_2", 5.0}};
  const std::vector<unsigned> expectedElemsOnEachProc {4, 4};
  run_decomposition_two_procs("rcb", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, NonuniformWeightsTwoBlocksParmetis)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_two_blocks();
  const BlockWeightsMap weightForBlock {{"block_1", 1.0}, {"block_2", 2.0}};
  const std::vector<unsigned> expectedElemsOnEachProc {5, 3};
  //
  //   <-- block_1 --> <-- block_2 -->
  //  o---o---o---o---o---o---o---o---o
  //  | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 2 |
  //  o---o---o---o---o---o---o---o---o
  //                      ^
  //                      New decomp places split here, with weight of 6 on each proc
  //
  run_decomposition_two_procs("parmetis", weightForBlock, expectedElemsOnEachProc);
}

TEST_F(TestBlockWeights, NonuniformWeightsTwoBlocksRCB)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  set_up_1x1x8_mesh_two_blocks();
  const BlockWeightsMap weightForBlock {{"block_1", 1.0}, {"block_2", 2.0}};
  const std::vector<unsigned> expectedElemsOnEachProc {5, 3};
  //
  //   <-- block_1 --> <-- block_2 -->
  //  o---o---o---o---o---o---o---o---o
  //  | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 2 |
  //  o---o---o---o---o---o---o---o---o
  //                      ^
  //                      New decomp places split here, with weight of 6 on each proc
  //
  run_decomposition_two_procs("rcb", weightForBlock, expectedElemsOnEachProc);
}

}
