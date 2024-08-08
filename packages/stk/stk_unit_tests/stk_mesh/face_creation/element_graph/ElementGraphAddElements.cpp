#include <gtest/gtest.h>
#include <algorithm>
#include <map>
#include <vector>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "ElementGraphTester.hpp"

namespace
{

class ElemGraphAddElementsToEmptyGraphTester : public stk::unit_test_util::MeshTestFixture
{
protected:
  ElemGraphAddElementsToEmptyGraphTester() : elementGraph(nullptr) { }
  ~ElemGraphAddElementsToEmptyGraphTester()
  {
    delete elementGraph;
  }

  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    create_elem_graph();
    stk::io::fill_mesh("generated:1x1x4", get_bulk());
    expect_graph_correct_after_moving_everything_to_proc0();
  }

  void expect_graph_correct_after_moving_everything_to_proc0()
  {
    add_elements_to_graph();
    test_edges();
  }

  void create_elem_graph()
  {
    elementGraph = new ElemElemGraphTester(get_bulk());
    EXPECT_EQ(0u, elementGraph->size());
  }

  void add_elements_to_graph()
  {
    stk::mesh::EntityVector elements_to_add;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEMENT_RANK, get_bulk().mesh_meta_data().locally_owned_part(), elements_to_add);
    elementGraph->add_elements(elements_to_add);
    EXPECT_EQ(elements_to_add.size(), elementGraph->size());
  }

  void test_num_edges()
  {
    std::vector<size_t> numGraphEdgesPerProc = {1, 2, 2, 1};
    size_t numGraphEdges = numGraphEdgesPerProc[get_bulk().parallel_rank()];
    EXPECT_EQ(numGraphEdges, elementGraph->num_edges());
    EXPECT_EQ(numGraphEdges, elementGraph->num_parallel_edges());
  }

  void test_edges()
  {
    test_num_edges();

    GraphEdges graphEdgesThisProc = goldElement1ToElement2SideOrdinalsPerOwningProc[get_bulk().parallel_rank()];
    for(GraphEdgeMock const &graphEdge : graphEdgesThisProc)
      EXPECT_EQ(graphEdge.sideOrdinalConnectingElement1ToElement2, elementGraph->get_side_from_element1_to_element2(graphEdge.element1, graphEdge.element2));
  }
protected:
  ElemElemGraphTester *elementGraph;

  std::vector<GraphEdges> goldElement1ToElement2SideOrdinalsPerOwningProc =
  {
    {{1, 2, 5}, {1, 3, -1}, {1, 4, -1}},
    {{2, 1, 4}, {2, 3,  5}, {2, 4, -1}},
    {{3, 2, 4}, {3, 4,  5}, {3, 1, -1}},
    {{4, 3, 4}, {4, 1, -1}, {4, 2, -1}}
  };
};
TEST_F(ElemGraphAddElementsToEmptyGraphTester, withAura)
{
  run_test_on_num_procs(4, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphAddElementsToEmptyGraphTester, withoutAura)
{
  run_test_on_num_procs(4, stk::mesh::BulkData::NO_AUTO_AURA);
}

struct NodeSharingInfo
{
  unsigned owningProc;
  unsigned sharedNodeId;
  unsigned procSharedTo;
};

void setup_node_sharing(stk::mesh::BulkData &mesh, const std::vector<NodeSharingInfo> & sharedNodes);

class ElemGraphAddElementsToExistingGraphTester : public stk::unit_test_util::MeshTestFixture
{
protected:
  ElemGraphAddElementsToExistingGraphTester() : elementGraph(nullptr), hexPart(nullptr) { }
  ~ElemGraphAddElementsToExistingGraphTester()
  {
    delete elementGraph;
  }

  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    hexPart = &get_meta().declare_part_with_topology("hex_part", stk::topology::HEX_8);
    stk::io::fill_mesh("generated:1x1x2", get_bulk());
    create_element_graph();
    expect_graph_correct_after_adding_two_new_elements();
  }

  void create_element_graph()
  {
    elementGraph = new ElemElemGraphTester(get_bulk());
    EXPECT_EQ(1u, elementGraph->size());
  }

  void expect_graph_correct_after_adding_two_new_elements()
  {
    declare_new_hex_element_per_proc();
    add_new_hex_element_to_existing_graph();
    test_edges();
  }

  void declare_new_hex_element_per_proc()
  {
    get_bulk().modification_begin();
    declare_hex_element();
    setup_node_sharing_for_added_element();
    get_bulk().modification_end();
  }

  void declare_hex_element()
  {
    const int rank = get_bulk().parallel_rank();
    stk::mesh::EntityIdVector hexNodeIds[]
    {
      {  12, 11, 18, 17, 8, 7, 16, 15 },
      {   8,  7, 16, 15, 4, 3, 14, 13 }
    };
    stk::mesh::declare_element(get_bulk(), *hexPart, hexElemIdsByProc[rank], hexNodeIds[rank]);
  }

  void setup_node_sharing_for_added_element()
  {
    std::vector<NodeSharingInfo> sharedNodes
    {
      {0, 3, 1}, {0, 4, 1}, {0, 11, 1}, {0, 12, 1}, {0, 15, 1}, {0, 16, 1}, // proc 0
      {1, 3, 0}, {1, 4, 0}, {1, 11, 0}, {1, 12, 0}, {1, 15, 0}, {1, 16, 0}, // proc 1
    };
    setup_node_sharing(get_bulk(), sharedNodes );
  }

  void add_new_hex_element_to_existing_graph()
  {
    const int rank = get_bulk().parallel_rank();
    stk::mesh::Entity elementToAdd = get_bulk().get_entity(stk::topology::ELEM_RANK, hexElemIdsByProc[rank]);
    elementGraph->add_elements({elementToAdd});
    EXPECT_EQ(2u, elementGraph->size());
  }

  void test_num_edges()
  {
    EXPECT_EQ(4u, elementGraph->num_edges());
    EXPECT_EQ(4u, elementGraph->num_parallel_edges());
  }

  void test_edges()
  {
    test_num_edges();

    GraphEdges graphEdgesThisProc = goldElement1ToElement2SideOrdinalsPerOwningProc[get_bulk().parallel_rank()];
    for(GraphEdgeMock const &graphEdge : graphEdgesThisProc)
      EXPECT_EQ(graphEdge.sideOrdinalConnectingElement1ToElement2, elementGraph->get_side_from_element1_to_element2(graphEdge.element1, graphEdge.element2));
  }

protected:
  ElemElemGraphTester *elementGraph;
  stk::mesh::Part *hexPart;

  stk::mesh::EntityIdVector hexElemIdsByProc = {3, 4};
  std::vector<GraphEdges> goldElement1ToElement2SideOrdinalsPerOwningProc =
  {
    {{1, 2, 5}, {1, 4, 2}, {3, 2, 0}, {3, 4, 5}},
    {{2, 1, 4}, {2, 3, 2}, {4, 1, 0}, {4, 3, 4}},
  };
};
TEST_F(ElemGraphAddElementsToExistingGraphTester, aura_on)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphAddElementsToExistingGraphTester, aura_off)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}

// THIS CODE EXISTS IN UnitTestSkinMesh.cpp AS WELL...doh!
void setup_node_sharing(stk::mesh::BulkData &mesh, const std::vector<NodeSharingInfo> & sharedNodes)
{
  for(const NodeSharingInfo &nodeSharing : sharedNodes)
  {
    if (static_cast<unsigned>(mesh.parallel_rank()) == nodeSharing.owningProc)
    {
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeSharing.sharedNodeId);
      mesh.add_node_sharing(node, nodeSharing.procSharedTo);
    }
  }
}

}
