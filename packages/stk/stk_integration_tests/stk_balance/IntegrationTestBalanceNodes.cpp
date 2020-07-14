#include "gtest/gtest.h"
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

class BalanceNodes : public stk::unit_test_util::MeshFixture {};

TEST_F(BalanceNodes, twoHex_initiallyImbalanced)
{
  if (get_parallel_size() != 2) return;

  stk::balance::StkBalanceSettings balanceSettings;
  balanceSettings.setUseNodeBalancer(true);

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::unit_test_util::setup_text_mesh(get_bulk(),
                                       "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                                       "0,2,HEX_8,5,6,7,8,9,10,11,12");

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());
  stk::balance::balanceStkMeshNodes(balanceSettings, get_bulk());

  const size_t numNodes = count_selected_entities(get_meta().locally_owned_part(),
                                                  get_bulk().buckets(stk::topology::NODE_RANK));

  EXPECT_EQ(numNodes, 6u);
}

TEST_F(BalanceNodes, twoHex_initiallyBalanced)
{
  if (get_parallel_size() != 2) return;

  stk::balance::StkBalanceSettings balanceSettings;
  balanceSettings.setUseNodeBalancer(true);

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::unit_test_util::setup_text_mesh(get_bulk(),
                                       "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                                       "1,2,HEX_8,5,6,7,8,9,10,11,12");

  stk::balance::balanceStkMeshNodes(balanceSettings, get_bulk());

  const size_t numNodes = count_selected_entities(get_meta().locally_owned_part(),
                                                  get_bulk().buckets(stk::topology::NODE_RANK));

  EXPECT_EQ(numNodes, 6u);
}

TEST_F(BalanceNodes, threeHex)
{
  if (get_parallel_size() != 2) return;

  stk::balance::StkBalanceSettings balanceSettings;
  balanceSettings.setUseNodeBalancer(true);

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::unit_test_util::setup_text_mesh(get_bulk(),
                                       "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                                       "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                                       "0,3,HEX_8,9,10,11,12,13,14,15,16");

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());
  stk::balance::balanceStkMeshNodes(balanceSettings, get_bulk());

  const size_t numNodes = count_selected_entities(get_meta().locally_owned_part(),
                                                  get_bulk().buckets(stk::topology::NODE_RANK));

  EXPECT_EQ(numNodes, 8u);
}

double get_node_imbalance(const stk::mesh::BulkData & bulk)
{
  stk::mesh::Selector localSelector = bulk.mesh_meta_data().locally_owned_part();
  stk::mesh::EntityVector ownedNodes;
  bulk.get_entities(stk::topology::NODE_RANK, localSelector, ownedNodes);

  const size_t numLocallyOwnedNodes = ownedNodes.size();
  size_t maxLocallyOwned = 0;
  size_t minLocallyOwned = 0;
  stk::all_reduce_max(bulk.parallel(), &numLocallyOwnedNodes, &maxLocallyOwned, 1);
  stk::all_reduce_min(bulk.parallel(), &numLocallyOwnedNodes, &minLocallyOwned, 1);
  return double(maxLocallyOwned) / double(minLocallyOwned);
}

TEST_F(BalanceNodes, bigMesh)
{
  if (get_parallel_size() != 4) return;

  setup_mesh("generated:16x16x16", stk::mesh::BulkData::NO_AUTO_AURA);

  const double targetLoadBalance = 1.01;
  stk::balance::StkBalanceSettings balanceSettings;
  balanceSettings.setUseNodeBalancer(true);
  balanceSettings.setNodeBalancerTargetLoadBalance(targetLoadBalance);
  balanceSettings.setNodeBalancerMaxIterations(10);
  stk::balance::balanceStkMeshNodes(balanceSettings, get_bulk());

  EXPECT_LT(get_node_imbalance(get_bulk()), targetLoadBalance);
}

