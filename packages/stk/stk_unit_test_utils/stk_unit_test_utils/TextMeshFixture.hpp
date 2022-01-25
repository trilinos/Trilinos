#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTTEXTMESH_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTTEXTMESH_HPP_

#include "TextMeshStkTopologyMapping.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field, etc
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <string>
#include <vector>

using EntityId = stk::mesh::EntityId;
using EntityIdVector = std::vector<EntityId>;
using Topology = stk::topology;

namespace stk
{
namespace unit_test_util
{
class TextMeshFixture : public stk::unit_test_util::MeshFixture
{
 protected:
  TextMeshFixture(unsigned spatialDim);

  void setup_text_mesh(const std::string& meshDesc);

  void setup_text_mesh(const std::string& meshDesc, const std::vector<double>& coordinates);

  void verify_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int sharingProc);

  void verify_num_elements(size_t goldCount);

  void verify_single_element(
      stk::mesh::EntityId elemId, const std::string& textMeshTopologyName, const stk::mesh::EntityIdVector& nodeIds);

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

}  // namespace unit_test_util
}  // namespace stk

#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTTEXTMESH_HPP_ */
