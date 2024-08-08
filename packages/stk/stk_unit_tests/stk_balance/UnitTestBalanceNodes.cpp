#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_balance/balance.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace {

class TestBalanceNodes : public stk::unit_test_util::MeshFixture
{
protected:
};

TEST_F(TestBalanceNodes, withEdges)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numProcs < 2 || numProcs > 4) { GTEST_SKIP(); }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::Part& edgePart = get_meta().declare_part("edges", stk::topology::EDGE_RANK);

  std::string meshDesc("generated:8x8x8");
  stk::io::fill_mesh(meshDesc, get_bulk());
  const int totalNumNodes = 9*9*9;

  stk::mesh::create_edges(get_bulk(), get_meta().universal_part(), &edgePart);

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rib");
  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  balanceSettings.setUseNodeBalancer(true);
  balanceSettings.setNodeBalancerTargetLoadBalance(1);
  balanceSettings.setNodeBalancerMaxIterations(4);
  stk::balance::balanceStkMeshNodes(balanceSettings, get_bulk());

  const int nodeCount = stk::mesh::count_entities(get_bulk(), stk::topology::NODE_RANK, get_meta().locally_owned_part());

  const int expectedNumOwnedNodesPerProc = totalNumNodes / numProcs;
  EXPECT_LE(std::abs(nodeCount - expectedNumOwnedNodesPerProc), 1);  // Within one of balanced
}

TEST_F(TestBalanceNodes, emptyProcessor)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numProcs < 2 || numProcs > 4) { GTEST_SKIP(); }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  std::string meshDesc("textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8");
  stk::io::fill_mesh(meshDesc, get_bulk());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setUseNodeBalancer(true);
  balanceSettings.setNodeBalancerTargetLoadBalance(1);
  balanceSettings.setNodeBalancerMaxIterations(4);

  stk::balance::balanceStkMeshNodes(balanceSettings, get_bulk());

  const int nodeCount = stk::mesh::count_entities(get_bulk(), stk::topology::NODE_RANK, get_meta().locally_owned_part());

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    EXPECT_EQ(nodeCount, 8);
  }
  else {
    EXPECT_EQ(nodeCount, 0);
  }
}

} // namespace
