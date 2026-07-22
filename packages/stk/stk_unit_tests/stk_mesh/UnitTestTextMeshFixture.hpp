#ifndef STK_STK_UNIT_TESTS_STK_MESH_UNITTESTTEXTMESHFIXTURE_HPP_
#define STK_STK_UNIT_TESTS_STK_MESH_UNITTESTTEXTMESHFIXTURE_HPP_

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>  // for Field
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field, etc
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/TextMeshFixture.hpp>
#include <stk_unit_test_utils/TextMeshStkTopologyMapping.hpp>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "mpi.h"

namespace
{
class TestTextMesh : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh() : TextMeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

class TestTextMeshAura : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMeshAura() : TextMeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

class TestTextMesh2d : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh2d() : TextMeshFixture(2)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

class TestTextMeshAura2d : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMeshAura2d() : TextMeshFixture(2)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

class TestTextMesh1d : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh1d() : TextMeshFixture(1)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

class TestTextMeshGraph : public stk::unit_test_util::TextMeshFixture
{
 protected:
  TestTextMeshGraph() : TextMeshFixture(3) {}

  class TextMeshGraph : public SideAdjacencyGraph
  {
  public:
    TextMeshGraph(const TextMeshData& data) : m_data(data) {}

    size_t get_num_elements() const override
    {
      return m_data.elementDataVec.size();
    }

    int get_element_proc(const size_t elemIndex) const override
    {
      const ElementData &elemData = m_data.elementDataVec[elemIndex];
      return elemData.proc;
    }

    bool element_has_any_node_on_proc(const size_t elemIndex, int proc) const override
    {
      const ElementData &elemData = m_data.elementDataVec[elemIndex];

      for (const EntityId &nodeId : elemData.nodeIds) {
        const std::set<int> &procsForNode = m_data.procs_for_node(nodeId);
        if (procsForNode.count(proc) > 0) {
          return true;
        }
      }

      return false;
    }

    const std::string& get_element_block_name(const size_t elemIndex) const override
    {
      const ElementData &elemData = m_data.elementDataVec[elemIndex];
      return elemData.partName;
    }

    const std::vector<EntityId>& get_element_node_ids(const size_t elemIndex) const override
    {
      const ElementData &elemData = m_data.elementDataVec[elemIndex];
      return elemData.nodeIds;
    }

    const Topology& get_element_topology(const size_t elemIndex) const override
    {
      const ElementData &elemData = m_data.elementDataVec[elemIndex];
      return elemData.topology;
    }

    EntityId get_element_id(const size_t elemIndex) const override
    {
      const ElementData &elemData = m_data.elementDataVec[elemIndex];
      return elemData.identifier;
    }

  private:
    const TextMeshData& m_data;
  };

  void dump_graph(std::ostream& out = std::cout) { m_graph->dump(m_data.elementDataVec, out); }

  void setup_text_mesh_graph(const std::string& meshDesc,
      const std::vector<std::string>& selectedBlocks = {},
      int proc = SideAdjacencyGraph::ANY_PROC)
  {
    TextMeshParser parser;
    m_data = parser.parse(meshDesc);
    m_graph = std::make_shared<TextMeshGraph>(m_data);
    m_graph->create_graph(selectedBlocks, proc);
  }

  void verify_side_adjacency(const std::vector<Adjacency>& goldNeighbors)
  {
    EXPECT_EQ(m_graph->size(), goldNeighbors.size());
    for (size_t i = 0; i < goldNeighbors.size(); ++i) {
      const auto& graphNeighborIndices = (*m_graph)[goldNeighbors[i].elemIndex];
      const auto& goldNeighborIndices = goldNeighbors[i].neighborIndices;

      unsigned numActualGoldConnections = 0;
      for (const auto& entry : goldNeighborIndices) {
        if (entry.second >= 0) {
          numActualGoldConnections++;
        }
      }

      EXPECT_EQ(numActualGoldConnections, graphNeighborIndices.connections.size());

      for (const auto& entry : goldNeighborIndices) {
        int side = entry.first + 1;
        SideAdjacencyGraph::IndexType neighborElemIndex = entry.second;

        if (neighborElemIndex >= 0) {
          EXPECT_LT(0, graphNeighborIndices.sideReference[side - 1]);
          EXPECT_TRUE(graphNeighborIndices.has_any_connection(side, neighborElemIndex));
        } else {
          EXPECT_EQ(0, graphNeighborIndices.sideReference[side - 1]);
          EXPECT_FALSE(graphNeighborIndices.has_any_connection(side));
        }
      }
    }
  }

  TextMeshData m_data;
  std::shared_ptr<TextMeshGraph> m_graph;
};

using TestTextMeshSkin = TestTextMesh;
}



#endif /* STK_STK_UNIT_TESTS_STK_MESH_UNITTESTTEXTMESHFIXTURE_HPP_ */
