#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTTEXTMESH_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTTEXTMESH_HPP_

#include "TextMeshStkTopologyMapping.hpp"
#include "TextMeshSideset.hpp"
#include "TextMeshNodeset.hpp"
#include "TextMeshUtils.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field, etc
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <string>
#include <vector>

using Topology = stk::unit_test_util::StkTopologyMapEntry;
using TopologyMapping = stk::unit_test_util::StkTopologyMapping;
using EntityId = stk::mesh::EntityId;
using EntityIdVector = std::vector<EntityId>;
using TextMeshData = text_mesh::TextMeshData<EntityId, Topology>;
using ElementData = text_mesh::ElementData<EntityId, Topology>;
using SidesetData = text_mesh::SidesetData<EntityId, Topology>;
using NodesetData = text_mesh::NodesetData<EntityId>;
using Coordinates = text_mesh::Coordinates<EntityId>;
using TextMeshParser = text_mesh::TextMeshParser<EntityId, TopologyMapping>;
using SideAdjacencyGraph = text_mesh::SideAdjacencyGraph<EntityId, Topology>;
using SideBlockInfo = text_mesh::SideBlockInfo;
using SideEntry = std::pair<EntityId, int>;
using SideVector = std::vector<SideEntry>;
using SplitType = text_mesh::SplitType;

struct Adjacency {
  using NeighborVector = std::vector<std::pair<int, SideAdjacencyGraph::IndexType>>;
  using SimpleNeighborVector = std::vector<SideAdjacencyGraph::IndexType>;

  size_t elemIndex;
  NeighborVector neighborIndices;

  Adjacency(size_t elemIndex_, const NeighborVector& neighborIndices_)
      : elemIndex(elemIndex_), neighborIndices(neighborIndices_)
  {
  }

  Adjacency(size_t elemIndex_, const SimpleNeighborVector& neighborIndices_)
      : elemIndex(elemIndex_), neighborIndices(get_full_neighbor_vector(neighborIndices_))
  {
  }

  NeighborVector get_full_neighbor_vector(const SimpleNeighborVector& simpleNeighborVector)
  {
    NeighborVector fullNeighborVector;
    fullNeighborVector.reserve(simpleNeighborVector.size());

    for (unsigned i = 0; i < simpleNeighborVector.size(); i++) {
      fullNeighborVector.push_back(std::make_pair(static_cast<int>(i), simpleNeighborVector[i]));
    }
    return fullNeighborVector;
  }
};

namespace stk
{
namespace unit_test_util
{

class TextMeshFixture : public stk::unit_test_util::MeshFixture
{
 protected:
  TextMeshFixture(unsigned spatialDim);

  void setup_text_mesh(const std::string& meshDesc);

  void verify_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int sharingProc);

  void verify_num_elements(size_t goldCount);

  void verify_single_element(
      stk::mesh::EntityId elemId, const std::string& textMeshTopologyName, const stk::mesh::EntityIdVector& nodeIds);

  void verify_num_sidesets(size_t goldCount);

  void verify_single_sideset(const std::string& name, const unsigned id, const SideVector& elemSidePairs);

  void verify_single_sideset(const std::string& name,
      const unsigned id,
      const std::vector<std::string>& subsets,
      const SideVector& elemSidePairs);

  void verify_num_nodesets(size_t goldCount);

  void verify_single_nodeset(const std::string& name, const unsigned id, const EntityIdVector& nodeIds);

  void verify_num_assemblies(size_t goldCount);

  void verify_single_assembly(const std::string& name, const unsigned id, const std::vector<std::string>& goldMembers);

  struct PartInfo {
    std::string blockName;
    std::set<stk::mesh::EntityId> ids;
  };

  void verify_part_membership(const std::vector<PartInfo> golds);

  using PartNameId = std::pair<std::string, unsigned>;

  void verify_part_ids(const std::vector<PartNameId>& golds);

  void verify_coordinates(const stk::mesh::EntityIdVector& goldNodeIds, const std::vector<double>& goldCoordinates);

  std::string get_topology_name(const std::string& textMeshTopologyName);

 private:
  void verify_nodes_on_element(stk::mesh::Entity element, const stk::mesh::EntityIdVector& goldNodeIds);

  stk::mesh::EntityIdVector get_node_ids(const stk::mesh::EntityVector& nodes);

  void verify_part(stk::mesh::Part* blockPart);

  void verify_elements_on_part(stk::mesh::Part* blockPart, const std::set<stk::mesh::EntityId>& goldIds);

  void verify_sideset_subset(stk::mesh::Part* sidesetPart, const unsigned id, const std::vector<std::string>& subsets);

  bool element_has_side_on_ordinal(stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal ordinal);

  class CoordinateVerifier
  {
   public:
    CoordinateVerifier(
        const stk::mesh::BulkData& b, const stk::mesh::EntityIdVector& ids, const std::vector<double>& coords);

    void verify();

   private:
    void verify_num_nodes();

    const double* get_nodal_coordinates(const stk::mesh::EntityId& nodeId);

    const stk::mesh::Entity get_node(const stk::mesh::EntityId& nodeId);

    void verify_nodal_coordinates(
        const stk::mesh::EntityId& nodeId, const double* goldCoords, const double* nodalCoords);

    std::string error_message(const stk::mesh::EntityId& nodeId, unsigned coordIndex);

    const stk::mesh::BulkData& bulk;
    const stk::mesh::MetaData& meta;

    const unsigned spatialDim;

    const stk::mesh::EntityIdVector& goldNodeIds;
    const std::vector<double>& goldCoordinates;
  };

  StkTopologyMapping m_topologyMapping;
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
TextMeshFixture : public stk::unit_test_util::MeshFixture
{
 protected:
  TextMeshFixture(unsigned spatialDim);

  void setup_text_mesh(const std::string& meshDesc);

  void verify_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int sharingProc);

  void verify_num_elements(size_t goldCount);

  void verify_single_element(
      stk::mesh::EntityId elemId, const std::string& textMeshTopologyName, const stk::mesh::EntityIdVector& nodeIds);

  void verify_num_sidesets(size_t goldCount);

  void verify_single_sideset(const std::string& name, const unsigned id, const SideVector& elemSidePairs);

  void verify_single_sideset(const std::string& name,
      const unsigned id,
      const std::vector<std::string>& subsets,
      const SideVector& elemSidePairs);

  void verify_num_nodesets(size_t goldCount);

  void verify_single_nodeset(const std::string& name, const unsigned id, const EntityIdVector& nodeIds);

  void verify_num_assemblies(size_t goldCount);

  void verify_single_assembly(const std::string& name, const unsigned id, const std::vector<std::string>& goldMembers);

  struct PartInfo {
    std::string blockName;
    std::set<stk::mesh::EntityId> ids;
  };

  void verify_part_membership(const std::vector<PartInfo> golds);

  using PartNameId = std::pair<std::string, unsigned>;

  void verify_part_ids(const std::vector<PartNameId>& golds);

  void verify_coordinates(const stk::mesh::EntityIdVector& goldNodeIds, const std::vector<double>& goldCoordinates);

  std::string get_topology_name(const std::string& textMeshTopologyName);

 private:
  void verify_nodes_on_element(stk::mesh::Entity element, const stk::mesh::EntityIdVector& goldNodeIds);

  stk::mesh::EntityIdVector get_node_ids(const stk::mesh::EntityVector& nodes);

  void verify_part(stk::mesh::Part* blockPart);

  void verify_elements_on_part(stk::mesh::Part* blockPart, const std::set<stk::mesh::EntityId>& goldIds);

  void verify_sideset_subset(stk::mesh::Part* sidesetPart, const unsigned id, const std::vector<std::string>& subsets);

  bool element_has_side_on_ordinal(stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal ordinal);

  class CoordinateVerifier
  {
   public:
    CoordinateVerifier(
        const stk::mesh::BulkData& b, const stk::mesh::EntityIdVector& ids, const std::vector<double>& coords);

    void verify();

   private:
    void verify_num_nodes();

    const double* get_nodal_coordinates(const stk::mesh::EntityId& nodeId);

    const stk::mesh::Entity get_node(const stk::mesh::EntityId& nodeId);

    void verify_nodal_coordinates(
        const stk::mesh::EntityId& nodeId, const double* goldCoords, const double* nodalCoords);

    std::string error_message(const stk::mesh::EntityId& nodeId, unsigned coordIndex);

    const stk::mesh::BulkData& bulk;
    const stk::mesh::MetaData& meta;

    const unsigned spatialDim;

    const stk::mesh::EntityIdVector& goldNodeIds;
    const std::vector<double>& goldCoordinates;
  };

  stk::unit_test_util::StkTopologyMapping m_topologyMapping;
};

} // namespace simple_fields

}  // namespace unit_test_util
}  // namespace stk

#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTTEXTMESH_HPP_ */
