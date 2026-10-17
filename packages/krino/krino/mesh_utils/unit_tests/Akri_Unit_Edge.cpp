// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Edge.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <array>
#include <unordered_set>
#include <vector>

namespace krino
{

//--------------------------------------------------------------------------
// The Edge value type packs the two edge-node local offsets into a single
// 64-bit key. These pin the struct semantics that the hashing/comparison
// containers in the refinement code rely on -- no mesh required.
//--------------------------------------------------------------------------
TEST(Edge, defaultConstructed_isInvalid)
{
  Edge edge;
  EXPECT_FALSE(edge.is_valid());
  EXPECT_EQ(Edge::InvalidEdge, edge.value());
}

TEST(Edge, nonzeroValue_isValid)
{
  Edge edge(42u);
  EXPECT_TRUE(edge.is_valid());
  EXPECT_EQ(42u, edge.value());
}

TEST(Edge, equalityAndOrdering_followUnderlyingValue)
{
  Edge a(5u);
  Edge b(5u);
  Edge c(7u);

  EXPECT_TRUE(a == b);
  EXPECT_TRUE(a == 5u);
  EXPECT_TRUE(a != c);
  EXPECT_TRUE(a < c);
  EXPECT_FALSE(c < a);

  // The specialized std::hash just returns the value, so equal edges hash
  // equally and can be de-duplicated in an unordered_set.
  std::unordered_set<Edge> edges;
  edges.insert(a);
  edges.insert(b); // same value as a
  edges.insert(c);
  EXPECT_EQ(2u, edges.size());
}

//--------------------------------------------------------------------------
// Mesh-based edge encoding/decoding, built on a single regular tet.
//--------------------------------------------------------------------------
class RegularTetEdge : public StkMeshTetFixture
{
protected:
  RegularTetEdge()
  {
    RegularTet meshSpec;
    set_valid_proc_sizes_for_test({1});
    if (is_valid_proc_size_for_test())
      StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Entity get_element() { return get_owned_elements()[0]; }
};

// edge_from_edge_nodes encodes the two nodes (ordered by identifier) and
// get_edge_nodes / fill_edge_nodes decode them; the pair must round-trip.
TEST_F(RegularTetEdge, edgeFromNodes_roundTripsThroughGetEdgeNodes)
{
  if (!is_valid_proc_size_for_test())
    return;

  const stk::mesh::BulkData & mesh = get_mesh();
  const stk::mesh::Entity elem = get_element();
  const stk::mesh::Entity * const nodes = mesh.begin_nodes(elem);
  const stk::mesh::Entity node0 = nodes[0];
  const stk::mesh::Entity node1 = nodes[1];

  // The encoding is order-independent in its inputs (it sorts by identifier),
  // so the two argument orders must yield the same edge.
  const Edge edge = edge_from_edge_nodes(mesh, node0, node1);
  const Edge edgeReversed = edge_from_edge_nodes(mesh, node1, node0);
  EXPECT_TRUE(edge.is_valid());
  EXPECT_EQ(edge, edgeReversed);

  // Decoding recovers the same two node entities (the lower-id node first).
  const std::array<stk::mesh::Entity, 2> decoded = get_edge_nodes(edge);
  const stk::mesh::Entity lowerIdNode = (mesh.identifier(node0) < mesh.identifier(node1)) ? node0 : node1;
  const stk::mesh::Entity higherIdNode = (mesh.identifier(node0) < mesh.identifier(node1)) ? node1 : node0;
  EXPECT_EQ(lowerIdNode, decoded[0]);
  EXPECT_EQ(higherIdNode, decoded[1]);

  // fill_edge_nodes must agree with get_edge_nodes.
  std::vector<stk::mesh::Entity> filled;
  fill_edge_nodes(edge, filled);
  ASSERT_EQ(2u, filled.size());
  EXPECT_EQ(decoded[0], filled[0]);
  EXPECT_EQ(decoded[1], filled[1]);
}

// A tet has 6 topological edges; fill_entity_edges / get_edges_of_elements
// must produce that many unique, valid edges.
TEST_F(RegularTetEdge, entityAndElementEdges_areSixUniqueValidEdgesForTet)
{
  if (!is_valid_proc_size_for_test())
    return;

  const stk::mesh::BulkData & mesh = get_mesh();
  const stk::mesh::Entity elem = get_element();

  std::vector<Edge> entityEdges;
  fill_entity_edges(mesh, elem, entityEdges);
  ASSERT_EQ(6u, entityEdges.size());
  for (const Edge & edge : entityEdges)
    EXPECT_TRUE(edge.is_valid());

  // All six edges are distinct.
  const std::unordered_set<Edge> uniqueEdges(entityEdges.begin(), entityEdges.end());
  EXPECT_EQ(6u, uniqueEdges.size());

  // get_edges_of_elements sorts-and-uniques the same set, so the count matches.
  const std::vector<Edge> elementEdges = get_edges_of_elements(mesh, {elem});
  EXPECT_EQ(6u, elementEdges.size());
}

}
