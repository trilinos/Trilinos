#include <gtest/gtest.h>  // for AssertHelper, EXPECT_EQ, etc

#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <string>

#include "stk_tools/mesh_tools/DisjointSet.hpp"
namespace stk::experimental
{

TEST(DisjointSetTests, merge_trees_with_different_nodes_throws)
{
  using Entity = stk::mesh::Entity;
  NodeElemKey n5e1(Entity(5), Entity(1));
  NodeElemKey n6e2(Entity(6), Entity(2));
  DisjointSet djSet;
  djSet.insert(n5e1);
  djSet.insert(n6e2);
  EXPECT_ANY_THROW(djSet.merge_nodes(n5e1, n6e2));
}

TEST(DisjointSetTests, throws_with_identical_insert)
{
  using Entity = stk::mesh::Entity;
  NodeElemKey n5e1(Entity(5), Entity(1));
  DisjointSet djSet;
  djSet.insert(n5e1);
  EXPECT_ANY_THROW(djSet.insert(n5e1));
}

TEST(DisjointSetTests, find_throws_from_nonexistent_key)
{
  using Entity = stk::mesh::Entity;
  NodeElemKey n5e1(Entity(5), Entity(1));
  DisjointSet djSet;
  EXPECT_ANY_THROW(djSet.find_root(n5e1));
}

TEST(DisjointSetTests, merge_trees_disjoint_set)
{
  using Entity = stk::mesh::Entity;

  DisjointSet dsNodes;
  NodeElemKey n5e1(Entity(5), Entity(1));
  NodeElemKey n5e2(Entity(5), Entity(2));
  NodeElemKey n5e3(Entity(5), Entity(3));
  NodeElemKey n5e4(Entity(5), Entity(4));
  dsNodes.insert(n5e1);
  dsNodes.insert(n5e2);
  dsNodes.insert(n5e3);
  dsNodes.insert(n5e4);

  dsNodes.merge_nodes(n5e1, n5e4);
  dsNodes.merge_nodes(n5e4, n5e3);
  dsNodes.merge_nodes(n5e3, n5e2);

  EXPECT_EQ(dsNodes.find_root(n5e1), dsNodes.find_root(n5e2));
  EXPECT_EQ(dsNodes.find_root(n5e2), dsNodes.find_root(n5e3));
  EXPECT_EQ(dsNodes.find_root(n5e3), dsNodes.find_root(n5e4));
}

struct DisjointSetFixture : public stk::unit_test_util::MeshFixture {
  DisjointSetFixture() : stk::unit_test_util::MeshFixture(3) { setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA); }
};

TEST_F(DisjointSetFixture, create_disjoint_set_and_merge_nodes_gen3x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  auto& meta = get_meta();
  auto& bulk = get_bulk();
  stk::io::fill_mesh("generated:3x2x2", bulk);
  DisjointSet disjointSet;

  auto owned = meta.universal_part() & meta.locally_owned_part();
  disjointSet.fill_set(bulk, owned);
  EXPECT_EQ(disjointSet.count_trees(), 8U * 12U);

  auto* currentRoot = &(disjointSet.begin()->second);
  for (auto& [key, node] : disjointSet) {
    if (currentRoot->get_node() == node.get_node()) {
      disjointSet.merge_nodes(currentRoot->key, key);
    } else {
      currentRoot = &node;
    }
  }
  EXPECT_EQ(disjointSet.count_trees(), 4U * 3U * 3U);
}

using KeyConnectivity = std::array<NodeElemKey, 8U>;
auto get_elem_connectivities(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& sel)
{
  std::vector<std::pair<stk::mesh::Entity, stk::mesh::ConnectedEntities>> origConns;
  std::vector<std::pair<stk::mesh::Entity, KeyConnectivity>> keyConns;
  auto elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, sel);
  for (const auto* bucket : elemBuckets) {
    for (const auto& elem : *bucket) {
      auto elemConn = bulk.get_connected_entities(elem, stk::topology::NODE_RANK);
      origConns.push_back(std::make_pair(elem, elemConn));
      KeyConnectivity keyConn_e;
      for (auto n = 0U; n < elemConn.size(); ++n) {
        keyConn_e[n].elem = elem;
        keyConn_e[n].node = elemConn[n];
      }
      keyConns.push_back(std::make_pair(elem, keyConn_e));
    }
  }
  return std::make_pair(origConns, keyConns);
}

TEST_F(DisjointSetFixture, fill_and_merge_reproduces_elem_conn_gen1x2x3)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  auto& meta = get_meta();
  auto& bulk = get_bulk();
  stk::io::fill_mesh("generated:1x2x3", bulk);
  DisjointSet disjointSet;

  auto owned = meta.universal_part() & meta.locally_owned_part();
  disjointSet.fill_set(bulk, owned);
  EXPECT_EQ(disjointSet.count_trees(), 8U * 6U);

  auto* currentRoot = &(disjointSet.begin()->second);
  for (auto& [key, node] : disjointSet) {
    if (currentRoot->get_node() == node.get_node()) {
      disjointSet.merge_nodes(currentRoot->key, key);
    } else {
      currentRoot = &node;
    }
  }
  EXPECT_EQ(disjointSet.count_trees(), 2U * 3U * 4U);

  auto [origConns, keyConns] = get_elem_connectivities(bulk, owned);
  ASSERT_EQ(origConns.size(), keyConns.size());
  for (auto e = 0U; e < origConns.size(); ++e) {
    const auto& [origElem, origConn_e] = origConns[e];
    const auto& [keyElem, keyConn_e] = keyConns[e];
    ASSERT_EQ(origElem, keyElem);
    ASSERT_EQ(origConn_e.size(), keyConn_e.size());
    for (auto n = 0U; n < origConn_e.size(); ++n) {
      auto newId = disjointSet.find_root(keyConn_e[n]).get_node_id();
      auto origId = bulk.identifier(origConn_e[n]);
      EXPECT_EQ(newId, origId) << "Node ID disagrees for element " << bulk.identifier(origElem) << ", node " << n;
    }
  }
}

}  // namespace stk::experimental
