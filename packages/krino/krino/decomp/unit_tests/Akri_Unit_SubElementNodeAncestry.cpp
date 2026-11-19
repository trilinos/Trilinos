// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <vector>

#include <gtest/gtest.h>

#include <Akri_CDMesh.hpp>
#include <Akri_Element.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_SubElementNodeAncestry.hpp>
#include <Akri_SubElementNodeOrder.hpp>
#include <Akri_Unit_DecompositionFixture.hpp>

namespace krino
{

namespace
{

class SubElementNodeAncestryFixture : public DecompositionFixture<RightTet, LSPerInterfacePolicy, 1>
{
public:
  SubElementNodeAncestryFixture()
  {
    set_valid_proc_sizes_for_test({1});
    build_mesh_from_spec(meshSpec);
    if (is_valid_proc_size_for_test())
    {
      setup_ls_fields();
      cdmesh->generate_nonconformal_elements();
    }
  }

  Mesh_Element & element() { return *cdmesh->elements.front(); }
  const SubElementNode * tet_node(const int i) { return element().get_nodes()[i]; }

  // The last of a chain of edge nodes, each on the edge from node 2 or node 3 of the tet
  // (alternately) to the one before, the first on the edge from node 0 to node 1, as cutting
  // crossings beside a sharp edge makes; the one before it, if wanted.
  const SubElementNode *
  nested_edge_node(const int depth, const SubElementNode ** previous = nullptr)
  {
    const SubElementNode * node =
        cdmesh->create_edge_node(&element(), tet_node(0), tet_node(1), 0.5);
    for (int level = 1; level < depth; ++level)
    {
      if (previous) *previous = node;
      node = cdmesh->create_edge_node(&element(), tet_node(2 + level % 2), node, 0.5);
    }
    return node;
  }
};

// Comparing ancestries walked every shared ancestor twice at each level, so its cost doubled with
// each level of nesting: at this depth, comparing a node with itself would not finish.
constexpr int depth = 60;

} // namespace

namespace
{
bool less_by_id(const SubElementNodeAncestry & a, const SubElementNodeAncestry & b)
{
  return less_by_global_id_because_it_doesnt_impact_structure(a, b);
}
} // namespace

// Mesh nodes come first, by id; child nodes are ordered by their parents, sorted.
TEST_F(SubElementNodeAncestryFixture, shallowNodes_orderedByAncestry)
{
  if (is_valid_proc_size_for_test())
  {
    const SubElementNode * child01 =
        cdmesh->create_edge_node(&element(), tet_node(0), tet_node(1), 0.5);
    const SubElementNode * child02 =
        cdmesh->create_edge_node(&element(), tet_node(2), tet_node(0), 0.5);
    const bool node0First = less_by_global_id_because_it_doesnt_impact_structure(*tet_node(0), *tet_node(1));
    EXPECT_EQ(node0First, less_by_id(tet_node(0)->get_ancestry(), tet_node(1)->get_ancestry()));
    EXPECT_TRUE(less_by_id(tet_node(0)->get_ancestry(), child01->get_ancestry()));
    EXPECT_FALSE(less_by_id(child01->get_ancestry(), tet_node(0)->get_ancestry()));
    const bool node1First = less_by_global_id_because_it_doesnt_impact_structure(*tet_node(1), *tet_node(2));
    EXPECT_EQ(node1First, less_by_id(child01->get_ancestry(), child02->get_ancestry()));
    EXPECT_EQ(!node1First, less_by_id(child02->get_ancestry(), child01->get_ancestry()));
  }
}

// Child nodes by their parents' ancestries, with the same order at every level: by coordinates,
// the child of nodes 0 and 2, at (0,0.5,0), is before the child of nodes 0 and 1, at (0.5,0,0),
// because node 2 (0,1,0) is before node 1 (1,0,0), whatever their ids.  The conditional order is by
// coordinates where the ids are not parallel consistent, and by id where they are.
TEST_F(SubElementNodeAncestryFixture, childNodes_orderedByTheirParentsWithTheSameOrderAtEveryLevel)
{
  if (is_valid_proc_size_for_test())
  {
    const SubElementNode * child01 =
        cdmesh->create_edge_node(&element(), tet_node(0), tet_node(1), 0.5);
    const SubElementNode * child02 =
        cdmesh->create_edge_node(&element(), tet_node(2), tet_node(0), 0.5);
    ASSERT_TRUE(less_by_global_id_because_it_doesnt_impact_structure(*tet_node(1), *tet_node(2)));

    EXPECT_TRUE(less_by_coordinates_because_it_impacts_structure(child02->get_ancestry(), child01->get_ancestry()));
    EXPECT_FALSE(less_by_coordinates_because_it_impacts_structure(child01->get_ancestry(), child02->get_ancestry()));
    EXPECT_TRUE(less_by_id(child01->get_ancestry(), child02->get_ancestry()));

    constexpr bool idsAreParallelConsistent = true;
    EXPECT_TRUE(less_conditionally_by_coordinates_because_it_impacts_structure(
        !idsAreParallelConsistent, child02->get_ancestry(), child01->get_ancestry()));
    EXPECT_FALSE(less_conditionally_by_coordinates_because_it_impacts_structure(
        !idsAreParallelConsistent, child01->get_ancestry(), child02->get_ancestry()));
    EXPECT_TRUE(less_conditionally_by_coordinates_because_it_impacts_structure(
        idsAreParallelConsistent, child01->get_ancestry(), child02->get_ancestry()));
  }
}

TEST_F(SubElementNodeAncestryFixture, deeplyNestedNodeComparedWithItself_notLess)
{
  if (is_valid_proc_size_for_test())
  {
    const SubElementNode * node = nested_edge_node(depth);
    EXPECT_FALSE(less_by_id(node->get_ancestry(), node->get_ancestry()));
    EXPECT_FALSE(less_by_coordinates_because_it_impacts_structure(node->get_ancestry(), node->get_ancestry()));
  }
}

TEST_F(SubElementNodeAncestryFixture, distinctDeeplyNestedNodes_strictlyOrdered)
{
  if (is_valid_proc_size_for_test())
  {
    const SubElementNode * previous = nullptr;
    const SubElementNode * node = nested_edge_node(depth, &previous);
    const bool nodeIsLess = less_by_id(node->get_ancestry(), previous->get_ancestry());
    const bool previousIsLess = less_by_id(previous->get_ancestry(), node->get_ancestry());
    EXPECT_NE(nodeIsLess, previousIsLess);
  }
}

// Making a subelement orders its nodes by their ancestries (to permute them consistently), which
// compares each node with itself.
TEST_F(SubElementNodeAncestryFixture, subelementWithDeeplyNestedNodes_made)
{
  if (is_valid_proc_size_for_test())
  {
    const SubElementNode * previous = nullptr;
    const SubElementNode * node = nested_edge_node(depth, &previous);
    const SubElement_Tet_4 subelement(
        {tet_node(0), tet_node(1), previous, node}, {0, 1, 2, 3}, &element());
    EXPECT_EQ(4u, subelement.get_nodes().size());
  }
}

} // namespace krino
