#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>  // for Field
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field, etc
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/TextMeshFixture.hpp>
#include <string>
#include <vector>
#include <strings.h>

#include "gtest/gtest.h"
#include "mpi.h"

namespace stk
{
namespace unit_test_util
{
namespace {
struct PartLessByName {
  inline bool operator()( const mesh::Part * lhs , const mesh::Part * rhs ) const
    { return strcasecmp(lhs->name().c_str(), rhs->name().c_str()) < 0; }
};
}

TextMeshFixture::TextMeshFixture(unsigned spatialDim)
  : stk::unit_test_util::MeshFixture(spatialDim)
{
  m_topologyMapping.initialize_topology_map();
}

std::string TextMeshFixture::get_topology_name(const std::string& textMeshTopologyName)
{
  return m_topologyMapping.topology(textMeshTopologyName).name();
}

void TextMeshFixture::setup_text_mesh(const std::string& meshDesc)
{
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
}

void TextMeshFixture::verify_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int sharingProc)
{
  std::vector<size_t> counts;
  stk::mesh::count_entities(get_meta().globally_shared_part(), get_bulk(), counts);
  EXPECT_EQ(nodeIds.size(), counts[stk::topology::NODE_RANK]);

  for (stk::mesh::EntityId nodeId : nodeIds) {
    EXPECT_TRUE(get_bulk().in_shared(stk::mesh::EntityKey(stk::topology::NODE_RANK, nodeId), sharingProc));
  }
}

void TextMeshFixture::verify_num_elements(size_t goldCount)
{
  std::vector<size_t> counts;
  stk::mesh::count_entities(get_meta().universal_part(), get_bulk(), counts);
  EXPECT_EQ(goldCount, counts[stk::topology::ELEM_RANK]);
}

void TextMeshFixture::verify_single_element(
    stk::mesh::EntityId elemId, const std::string& textMeshTopologyName, const stk::mesh::EntityIdVector& nodeIds)
{
  stk::topology topology = m_topologyMapping.topology(textMeshTopologyName).topology;
  stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
  EXPECT_TRUE(get_bulk().is_valid(element));
  EXPECT_EQ(topology, get_bulk().bucket(element).topology());
  verify_nodes_on_element(element, nodeIds);
}

void TextMeshFixture::verify_num_sidesets(size_t goldCount)
{
  EXPECT_EQ(goldCount, get_bulk().get_number_of_sidesets());
}

void TextMeshFixture::verify_single_sideset(const std::string& name, const unsigned id, const SideVector& elemSidePairs)
{
  verify_single_sideset(name, id, std::vector<std::string>{}, elemSidePairs);
}

void TextMeshFixture::verify_sideset_subset(
    stk::mesh::Part* sidesetPart, const unsigned id, const std::vector<std::string>& subsetNames)
{
  EXPECT_TRUE(nullptr != sidesetPart);
  EXPECT_TRUE(stk::io::is_part_io_part(*sidesetPart));
  EXPECT_EQ(get_meta().side_rank(), sidesetPart->primary_entity_rank());
  EXPECT_EQ(id, sidesetPart->id());
  EXPECT_EQ(sidesetPart->subsets().size(), subsetNames.size());

  for (const std::string& subsetName : subsetNames) {
    stk::mesh::Part* subsetPart = get_meta().get_part(subsetName);
    EXPECT_TRUE(nullptr != subsetPart);
    EXPECT_TRUE(stk::io::is_part_io_part(*subsetPart));
    EXPECT_EQ(get_meta().side_rank(), subsetPart->primary_entity_rank());
    EXPECT_EQ(id, subsetPart->id());

    auto iter = std::find(sidesetPart->subsets().begin(), sidesetPart->subsets().end(), subsetPart);
    EXPECT_NE(iter, sidesetPart->subsets().end());
  }
}

bool TextMeshFixture::element_has_side_on_ordinal(stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal ordinal)
{
  stk::mesh::ConnectivityOrdinal const* ordinals = get_bulk().begin_ordinals(element, get_meta().side_rank());
  unsigned numOrdinals = get_bulk().num_connectivity(element, get_meta().side_rank());
  bool foundOrdinal = false;
  for (unsigned i = 0; i < numOrdinals; ++i) {
    if (ordinals[i] == ordinal) {
      foundOrdinal = true;
      break;
    }
  }
  return foundOrdinal;
}

void TextMeshFixture::verify_single_sideset(const std::string& name,
    const unsigned id,
    const std::vector<std::string>& subsets,
    const SideVector& elemSidePairs)
{
  stk::mesh::Part* sidesetPart = get_meta().get_part(name);
  verify_sideset_subset(sidesetPart, id, subsets);

  EXPECT_TRUE(get_bulk().does_sideset_exist(*sidesetPart));

  stk::mesh::SideSet& sideset = get_bulk().get_sideset(*sidesetPart);

  for (const SideEntry& sideEntry : elemSidePairs) {
    stk::mesh::EntityId elemId = sideEntry.first;
    int side = sideEntry.second;

    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    EXPECT_TRUE(get_bulk().is_valid(element));

    stk::topology topology = get_bulk().bucket(element).topology();
    EXPECT_TRUE(side > 0);
    EXPECT_TRUE(side <= (int) topology.num_sides());

    stk::mesh::ConnectivityOrdinal ordinal = side - 1;
    EXPECT_TRUE(sideset.contains(element, ordinal));

    EXPECT_TRUE(element_has_side_on_ordinal(element, ordinal));
  }
}

void TextMeshFixture::verify_num_nodesets(size_t goldCount)
{
  size_t numNodesets = 0;
  for (stk::mesh::Part* part : get_meta().get_parts()) {
    if (stk::io::is_part_io_part(*part) && part->primary_entity_rank() == stk::topology::NODE_RANK) {
      numNodesets++;
    }
  }
  EXPECT_EQ(goldCount, numNodesets);
}

void TextMeshFixture::verify_single_nodeset(
    const std::string& name, const unsigned id, const stk::mesh::EntityIdVector& nodeIds)
{
  stk::mesh::Part* nodesetPart = get_meta().get_part(name);
  EXPECT_TRUE(nullptr != nodesetPart);
  EXPECT_EQ(stk::topology::NODE_RANK, nodesetPart->primary_entity_rank());
  EXPECT_TRUE(stk::io::is_part_io_part(*nodesetPart));
  EXPECT_EQ(id, nodesetPart->id());

  for (stk::mesh::EntityId nodeId : nodeIds) {
    stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
    EXPECT_TRUE(get_bulk().is_valid(node));
    EXPECT_TRUE(get_bulk().bucket(node).member(*nodesetPart));
  }
}

void TextMeshFixture::verify_num_assemblies(size_t goldCount)
{
  EXPECT_EQ(goldCount, stk::io::get_assembly_names(get_meta()).size());
}

void TextMeshFixture::verify_single_assembly(const std::string& name, const unsigned id,
                                             const std::vector<std::string>& goldMembers)
{
  stk::mesh::Part* assemblyPart = get_meta().get_part(name);
  EXPECT_TRUE(nullptr != assemblyPart);
  EXPECT_TRUE(stk::io::is_part_assembly_io_part(*assemblyPart));
  EXPECT_EQ(id, assemblyPart->id());

  stk::mesh::PartVector leafParts = stk::io::get_unique_leaf_parts(get_meta(), name);
  std::sort(leafParts.begin(), leafParts.end(), PartLessByName());
  EXPECT_EQ(goldMembers.size(), leafParts.size());

  for(size_t i=0; i<goldMembers.size(); i++) {
    const std::string& goldMember = goldMembers[i];
    stk::mesh::Part* goldPart = get_meta().get_part(goldMember);
    EXPECT_TRUE(nullptr != goldPart);
    EXPECT_EQ(goldPart->mesh_meta_data_ordinal(), leafParts[i]->mesh_meta_data_ordinal())
      << "Comparison failure for " << name << ": " << goldPart->name() << " <-> " << leafParts[i]->name();
  }

}

void TextMeshFixture::verify_part_membership(const std::vector<PartInfo> golds)
{
  for (const PartInfo& gold : golds) {
    stk::mesh::Part* blockPart = get_meta().get_part(gold.blockName);

    verify_part(blockPart);
    verify_elements_on_part(blockPart, gold.ids);
  }
}

void TextMeshFixture::verify_part_ids(const std::vector<PartNameId>& golds)
{
  for (const PartNameId& gold : golds) {
    stk::mesh::Part* blockPart = get_meta().get_part(gold.first);

    verify_part(blockPart);
    EXPECT_EQ(blockPart->id(), gold.second);
  }
}

void TextMeshFixture::verify_coordinates(
    const stk::mesh::EntityIdVector& goldNodeIds, const std::vector<double>& goldCoordinates)
{
  CoordinateVerifier cv(get_bulk(), goldNodeIds, goldCoordinates);
  cv.verify();
}

void TextMeshFixture::verify_nodes_on_element(stk::mesh::Entity element, const stk::mesh::EntityIdVector& goldNodeIds)
{
  stk::mesh::EntityVector nodes(get_bulk().begin_nodes(element), get_bulk().end_nodes(element));
  EXPECT_EQ(goldNodeIds, get_node_ids(nodes));
}

stk::mesh::EntityIdVector TextMeshFixture::get_node_ids(const stk::mesh::EntityVector& nodes)
{
  stk::mesh::EntityIdVector nodeIds;
  for (const stk::mesh::Entity& node : nodes) {
    nodeIds.emplace_back(get_bulk().identifier(node));
  }
  return nodeIds;
}

void TextMeshFixture::verify_part(stk::mesh::Part* blockPart)
{
  ASSERT_TRUE(blockPart != nullptr);
  EXPECT_TRUE(stk::io::is_part_io_part(*blockPart));
}

void TextMeshFixture::verify_elements_on_part(stk::mesh::Part* blockPart, const std::set<stk::mesh::EntityId>& goldIds)
{
  stk::mesh::EntityVector elems;
  stk::mesh::get_selected_entities(*blockPart, get_bulk().buckets(stk::topology::ELEM_RANK), elems);

  ASSERT_EQ(goldIds.size(), elems.size());
  for (const stk::mesh::Entity& elem : elems) {
    stk::mesh::EntityId elemId = get_bulk().identifier(elem);
    EXPECT_EQ(1u, goldIds.count(elemId));
  }
}

TextMeshFixture::CoordinateVerifier::CoordinateVerifier(const stk::mesh::BulkData& b,
                                                        const stk::mesh::EntityIdVector& ids,
                                                        const std::vector<double>& coords)
  : bulk(b),
    meta(bulk.mesh_meta_data()),
    spatialDim(meta.spatial_dimension()),
    goldNodeIds(ids),
    goldCoordinates(coords)
{
}

void TextMeshFixture::CoordinateVerifier::verify()
{
  verify_num_nodes();

  for (size_t nodeIndex = 0; nodeIndex < goldNodeIds.size(); nodeIndex++) {
    stk::mesh::EntityId nodeId = goldNodeIds[nodeIndex];
    EXPECT_TRUE(bulk.is_valid(get_node(nodeId)));

    const double* nodalCoords = get_nodal_coordinates(nodeId);
    const double* goldCoords = &goldCoordinates[nodeIndex * spatialDim];

    verify_nodal_coordinates(nodeId, goldCoords, nodalCoords);
  }
}

void TextMeshFixture::CoordinateVerifier::verify_num_nodes()
{
  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, nodes);
  EXPECT_EQ(goldNodeIds.size(), nodes.size());
}

const double* TextMeshFixture::CoordinateVerifier::get_nodal_coordinates(const stk::mesh::EntityId& nodeId)
{
  const stk::mesh::Field<double>& coordsField =
      static_cast<const stk::mesh::Field<double>&>(*meta.coordinate_field());
  return stk::mesh::field_data(coordsField, get_node(nodeId));
}

const stk::mesh::Entity TextMeshFixture::CoordinateVerifier::get_node(const stk::mesh::EntityId& nodeId)
{
  return bulk.get_entity(stk::topology::NODE_RANK, nodeId);
}

void TextMeshFixture::CoordinateVerifier::verify_nodal_coordinates(
    const stk::mesh::EntityId& nodeId, const double* goldCoords, const double* nodalCoords)
{
  for (unsigned i = 0; i < spatialDim; i++) {
    EXPECT_NEAR(goldCoords[i], nodalCoords[i], 1.0e-9) << error_message(nodeId, i);
  }
}

std::string TextMeshFixture::CoordinateVerifier::error_message(const stk::mesh::EntityId& nodeId, unsigned coordIndex)
{
  std::stringstream message;
  message << "Proc " << bulk.parallel_rank() << ", Node ID " << nodeId << ", coord index " << coordIndex;
  return message.str();
}

namespace simple_fields {

TextMeshFixture::TextMeshFixture(unsigned spatialDim)
  : stk::unit_test_util::MeshFixture(spatialDim)
{
  m_topologyMapping.initialize_topology_map();
}

std::string TextMeshFixture::get_topology_name(const std::string& textMeshTopologyName)
{
  return m_topologyMapping.topology(textMeshTopologyName).name();
}

void TextMeshFixture::setup_text_mesh(const std::string& meshDesc)
{
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
}

void TextMeshFixture::verify_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int sharingProc)
{
  std::vector<size_t> counts;
  stk::mesh::count_entities(get_meta().globally_shared_part(), get_bulk(), counts);
  EXPECT_EQ(nodeIds.size(), counts[stk::topology::NODE_RANK]);

  for (stk::mesh::EntityId nodeId : nodeIds) {
    EXPECT_TRUE(get_bulk().in_shared(stk::mesh::EntityKey(stk::topology::NODE_RANK, nodeId), sharingProc));
  }
}

void TextMeshFixture::verify_num_elements(size_t goldCount)
{
  std::vector<size_t> counts;
  stk::mesh::count_entities(get_meta().universal_part(), get_bulk(), counts);
  EXPECT_EQ(goldCount, counts[stk::topology::ELEM_RANK]);
}

void TextMeshFixture::verify_single_element(
    stk::mesh::EntityId elemId, const std::string& textMeshTopologyName, const stk::mesh::EntityIdVector& nodeIds)
{
  stk::topology topology = m_topologyMapping.topology(textMeshTopologyName).topology;
  stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
  EXPECT_TRUE(get_bulk().is_valid(element));
  EXPECT_EQ(topology, get_bulk().bucket(element).topology());
  verify_nodes_on_element(element, nodeIds);
}

void TextMeshFixture::verify_num_sidesets(size_t goldCount)
{
  EXPECT_EQ(goldCount, get_bulk().get_number_of_sidesets());
}

void TextMeshFixture::verify_single_sideset(const std::string& name, const unsigned id, const SideVector& elemSidePairs)
{
  verify_single_sideset(name, id, std::vector<std::string>{}, elemSidePairs);
}

void TextMeshFixture::verify_sideset_subset(
    stk::mesh::Part* sidesetPart, const unsigned id, const std::vector<std::string>& subsetNames)
{
  EXPECT_TRUE(nullptr != sidesetPart);
  EXPECT_TRUE(stk::io::is_part_io_part(*sidesetPart));
  EXPECT_EQ(get_meta().side_rank(), sidesetPart->primary_entity_rank());
  EXPECT_EQ(id, sidesetPart->id());
  EXPECT_EQ(sidesetPart->subsets().size(), subsetNames.size());

  for (const std::string& subsetName : subsetNames) {
    stk::mesh::Part* subsetPart = get_meta().get_part(subsetName);
    EXPECT_TRUE(nullptr != subsetPart);
    EXPECT_TRUE(stk::io::is_part_io_part(*subsetPart));
    EXPECT_EQ(get_meta().side_rank(), subsetPart->primary_entity_rank());
    EXPECT_EQ(id, subsetPart->id());

    auto iter = std::find(sidesetPart->subsets().begin(), sidesetPart->subsets().end(), subsetPart);
    EXPECT_NE(iter, sidesetPart->subsets().end());
  }
}

bool TextMeshFixture::element_has_side_on_ordinal(stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal ordinal)
{
  stk::mesh::ConnectivityOrdinal const* ordinals = get_bulk().begin_ordinals(element, get_meta().side_rank());
  unsigned numOrdinals = get_bulk().num_connectivity(element, get_meta().side_rank());
  bool foundOrdinal = false;
  for (unsigned i = 0; i < numOrdinals; ++i) {
    if (ordinals[i] == ordinal) {
      foundOrdinal = true;
      break;
    }
  }
  return foundOrdinal;
}

void TextMeshFixture::verify_single_sideset(const std::string& name,
    const unsigned id,
    const std::vector<std::string>& subsets,
    const SideVector& elemSidePairs)
{
  stk::mesh::Part* sidesetPart = get_meta().get_part(name);
  verify_sideset_subset(sidesetPart, id, subsets);

  EXPECT_TRUE(get_bulk().does_sideset_exist(*sidesetPart));

  stk::mesh::SideSet& sideset = get_bulk().get_sideset(*sidesetPart);

  for (const SideEntry& sideEntry : elemSidePairs) {
    stk::mesh::EntityId elemId = sideEntry.first;
    int side = sideEntry.second;

    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    EXPECT_TRUE(get_bulk().is_valid(element));

    stk::topology topology = get_bulk().bucket(element).topology();
    EXPECT_TRUE(side > 0);
    EXPECT_TRUE(side <= (int) topology.num_sides());

    stk::mesh::ConnectivityOrdinal ordinal = side - 1;
    EXPECT_TRUE(sideset.contains(element, ordinal));

    EXPECT_TRUE(element_has_side_on_ordinal(element, ordinal));
  }
}

void TextMeshFixture::verify_num_nodesets(size_t goldCount)
{
  size_t numNodesets = 0;
  for (stk::mesh::Part* part : get_meta().get_parts()) {
    if (stk::io::is_part_io_part(*part) && part->primary_entity_rank() == stk::topology::NODE_RANK) {
      numNodesets++;
    }
  }
  EXPECT_EQ(goldCount, numNodesets);
}

void TextMeshFixture::verify_single_nodeset(
    const std::string& name, const unsigned id, const stk::mesh::EntityIdVector& nodeIds)
{
  stk::mesh::Part* nodesetPart = get_meta().get_part(name);
  EXPECT_TRUE(nullptr != nodesetPart);
  EXPECT_EQ(stk::topology::NODE_RANK, nodesetPart->primary_entity_rank());
  EXPECT_TRUE(stk::io::is_part_io_part(*nodesetPart));
  EXPECT_EQ(id, nodesetPart->id());

  for (stk::mesh::EntityId nodeId : nodeIds) {
    stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
    EXPECT_TRUE(get_bulk().is_valid(node));
    EXPECT_TRUE(get_bulk().bucket(node).member(*nodesetPart));
  }
}

void TextMeshFixture::verify_num_assemblies(size_t goldCount)
{
  EXPECT_EQ(goldCount, stk::io::get_assembly_names(get_meta()).size());
}

void TextMeshFixture::verify_single_assembly(const std::string& name, const unsigned id,
                                             const std::vector<std::string>& goldMembers)
{
  stk::mesh::Part* assemblyPart = get_meta().get_part(name);
  EXPECT_TRUE(nullptr != assemblyPart);
  EXPECT_TRUE(stk::io::is_part_assembly_io_part(*assemblyPart));
  EXPECT_EQ(id, assemblyPart->id());

  stk::mesh::PartVector leafParts = stk::io::get_unique_leaf_parts(get_meta(), name);
  std::sort(leafParts.begin(), leafParts.end(), PartLessByName());
  EXPECT_EQ(goldMembers.size(), leafParts.size());

  for(size_t i=0; i<goldMembers.size(); i++) {
    const std::string& goldMember = goldMembers[i];
    stk::mesh::Part* goldPart = get_meta().get_part(goldMember);
    EXPECT_TRUE(nullptr != goldPart);
    EXPECT_EQ(goldPart->mesh_meta_data_ordinal(), leafParts[i]->mesh_meta_data_ordinal())
      << "Comparison failure for " << name << ": " << goldPart->name() << " <-> " << leafParts[i]->name();
  }

}

void TextMeshFixture::verify_part_membership(const std::vector<PartInfo> golds)
{
  for (const PartInfo& gold : golds) {
    stk::mesh::Part* blockPart = get_meta().get_part(gold.blockName);

    verify_part(blockPart);
    verify_elements_on_part(blockPart, gold.ids);
  }
}

void TextMeshFixture::verify_part_ids(const std::vector<PartNameId>& golds)
{
  for (const PartNameId& gold : golds) {
    stk::mesh::Part* blockPart = get_meta().get_part(gold.first);

    verify_part(blockPart);
    EXPECT_EQ(blockPart->id(), gold.second);
  }
}

void TextMeshFixture::verify_coordinates(
    const stk::mesh::EntityIdVector& goldNodeIds, const std::vector<double>& goldCoordinates)
{
  CoordinateVerifier cv(get_bulk(), goldNodeIds, goldCoordinates);
  cv.verify();
}

void TextMeshFixture::verify_nodes_on_element(stk::mesh::Entity element, const stk::mesh::EntityIdVector& goldNodeIds)
{
  stk::mesh::EntityVector nodes(get_bulk().begin_nodes(element), get_bulk().end_nodes(element));
  EXPECT_EQ(goldNodeIds, get_node_ids(nodes));
}

stk::mesh::EntityIdVector TextMeshFixture::get_node_ids(const stk::mesh::EntityVector& nodes)
{
  stk::mesh::EntityIdVector nodeIds;
  for (const stk::mesh::Entity& node : nodes) {
    nodeIds.emplace_back(get_bulk().identifier(node));
  }
  return nodeIds;
}

void TextMeshFixture::verify_part(stk::mesh::Part* blockPart)
{
  ASSERT_TRUE(blockPart != nullptr);
  EXPECT_TRUE(stk::io::is_part_io_part(*blockPart));
}

void TextMeshFixture::verify_elements_on_part(stk::mesh::Part* blockPart, const std::set<stk::mesh::EntityId>& goldIds)
{
  stk::mesh::EntityVector elems;
  stk::mesh::get_selected_entities(*blockPart, get_bulk().buckets(stk::topology::ELEM_RANK), elems);

  ASSERT_EQ(goldIds.size(), elems.size());
  for (const stk::mesh::Entity& elem : elems) {
    stk::mesh::EntityId elemId = get_bulk().identifier(elem);
    EXPECT_EQ(1u, goldIds.count(elemId));
  }
}

TextMeshFixture::CoordinateVerifier::CoordinateVerifier(const stk::mesh::BulkData& b,
                                                        const stk::mesh::EntityIdVector& ids,
                                                        const std::vector<double>& coords)
  : bulk(b),
    meta(bulk.mesh_meta_data()),
    spatialDim(meta.spatial_dimension()),
    goldNodeIds(ids),
    goldCoordinates(coords)
{
}

void TextMeshFixture::CoordinateVerifier::verify()
{
  verify_num_nodes();

  for (size_t nodeIndex = 0; nodeIndex < goldNodeIds.size(); nodeIndex++) {
    stk::mesh::EntityId nodeId = goldNodeIds[nodeIndex];
    EXPECT_TRUE(bulk.is_valid(get_node(nodeId)));

    const double* nodalCoords = get_nodal_coordinates(nodeId);
    const double* goldCoords = &goldCoordinates[nodeIndex * spatialDim];

    verify_nodal_coordinates(nodeId, goldCoords, nodalCoords);
  }
}

void TextMeshFixture::CoordinateVerifier::verify_num_nodes()
{
  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, nodes);
  EXPECT_EQ(goldNodeIds.size(), nodes.size());
}

const double* TextMeshFixture::CoordinateVerifier::get_nodal_coordinates(const stk::mesh::EntityId& nodeId)
{
  const stk::mesh::Field<double>& coordsField =
      static_cast<const stk::mesh::Field<double>&>(*meta.coordinate_field());
  return stk::mesh::field_data(coordsField, get_node(nodeId));
}

const stk::mesh::Entity TextMeshFixture::CoordinateVerifier::get_node(const stk::mesh::EntityId& nodeId)
{
  return bulk.get_entity(stk::topology::NODE_RANK, nodeId);
}

void TextMeshFixture::CoordinateVerifier::verify_nodal_coordinates(
    const stk::mesh::EntityId& nodeId, const double* goldCoords, const double* nodalCoords)
{
  for (unsigned i = 0; i < spatialDim; i++) {
    EXPECT_NEAR(goldCoords[i], nodalCoords[i], 1.0e-9) << error_message(nodeId, i);
  }
}

std::string TextMeshFixture::CoordinateVerifier::error_message(const stk::mesh::EntityId& nodeId, unsigned coordIndex)
{
  std::stringstream message;
  message << "Proc " << bulk.parallel_rank() << ", Node ID " << nodeId << ", coord index " << coordIndex;
  return message.str();
}

} // namespace simple_fields

}  // namespace unit_test_util
}  // namespace stk
