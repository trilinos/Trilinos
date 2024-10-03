// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "ConstructedMesh.hpp"
#include <ctype.h>                                   // for toupper
#include <stddef.h>                                  // for size_t
#include <algorithm>                                 // for remove, etc
#include <iterator>                                  // for insert_iterator
#include <map>
#include <set>                                       // for set
#include <sstream>                                   // for operator<<, etc
#include <stk_io/IossBridge.hpp>                     // for is_part_io_part, etc
#include <stk_mesh/base/BulkData.hpp>                // for BulkData
#include <stk_mesh/base/FEMHelpers.hpp>              // for declare_element
#include <stk_mesh/base/Field.hpp>                   // for Field
#include <stk_mesh/base/GetEntities.hpp>             // for get_entities
#include <stk_mesh/base/MetaData.hpp>                // for MetaData, etc
#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector

#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"       // for Cartesian
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"               // for field_data
#include "stk_mesh/base/Types.hpp"                   // for EntityId, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/util/ReportHandler.hpp"           // for ThrowRequireMsg

namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{

void ConstructedMesh::create_block_elements_and_nodes(stk::mesh::BulkData& bulk,
                                                      const ConstructedElementBlock& block,
                                                      const unsigned elemIdOffset)
{
  stk::mesh::Part* part = bulk.mesh_meta_data().get_part(block.name);
  STK_ThrowRequire(nullptr != part);

  size_t elementIndex = elemIdOffset;

  for(size_t elementCount = 0; elementCount < block.connectivityIndex.size(); ++elementIndex, ++elementCount) {
    stk::mesh::EntityIdVector elementNodeIds;
    for(size_t i=0; i<block.connectivityIndex[elementCount].size(); ++i) {
      int nodeIndex = block.connectivityIndex[elementCount][i];
      ASSERT_TRUE(nodeIndex > 0);
      stk::mesh::EntityId nodeId = m_nodeIds[nodeIndex - 1];
      elementNodeIds.push_back(nodeId);
    }
    stk::mesh::declare_element(bulk, *part, m_elemIds[elementIndex], elementNodeIds);
  }
}

void ConstructedMesh::verify_mesh_data(const stk::mesh::MetaData& meta)
{
  ASSERT_TRUE(meta.spatial_dimension() == 2 || meta.spatial_dimension() == 3);
  ASSERT_EQ(m_spatialDimension, meta.spatial_dimension());

  ASSERT_EQ(m_xCoords.size(), m_nodeIds.size());
  ASSERT_EQ(m_yCoords.size(), m_nodeIds.size());
  if(m_spatialDimension == 3) {
    ASSERT_EQ(m_zCoords.size(), m_nodeIds.size());
  }
}

void ConstructedMesh::populate_bulk_data(stk::mesh::BulkData& bulk)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  for(const ConstructedElementBlock& elemBlock : m_elemBlocks) {
    stk::mesh::Part& block = meta.declare_part_with_topology(elemBlock.name, elemBlock.topology);
    stk::io::put_io_part_attribute(block);
    meta.set_part_id(block, elemBlock.id);
  }

  stk::mesh::Field<double> & coordsField = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates", 1);
  stk::mesh::put_field_on_mesh(coordsField, meta.universal_part(), m_spatialDimension, nullptr);
  stk::io::set_field_output_type(coordsField, stk::io::FieldOutputType::VECTOR_3D);

  bulk.modification_begin();
  if(bulk.parallel_rank() == 0) {
    unsigned elemIdOffset = 0;
    for(unsigned i=0; i<m_elemBlocks.size(); ++i) {
      for(const std::vector<unsigned>& connectivity : m_elemBlocks[i].connectivityIndex) {
        ASSERT_EQ(m_elemBlocks[i].topology.num_nodes(), connectivity.size());
      }
      create_block_elements_and_nodes(bulk, m_elemBlocks[i], elemIdOffset);
      elemIdOffset += m_elemBlocks[i].connectivityIndex.size();
    }

    for(size_t nodeIndex=0; nodeIndex < m_nodeIds.size(); nodeIndex++)
    {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, m_nodeIds[nodeIndex]);
      ASSERT_TRUE(bulk.is_valid(node)) << "Invalid node id: " << m_nodeIds[nodeIndex];
      double * nodalCoords = stk::mesh::field_data(coordsField, node);

      nodalCoords[0] = m_xCoords[nodeIndex];
      nodalCoords[1] = m_yCoords[nodeIndex];

      if(m_spatialDimension == 3) {
        nodalCoords[2] = m_zCoords[nodeIndex];
      }
    }
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec changeList;

  if(bulk.parallel_rank() == 0 && bulk.parallel_size() > 1) {
    for(const stk::mesh::EntityIdProc& entry : m_elemDistribution) {
      stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, entry.first);
      ASSERT_TRUE(bulk.is_valid(elem)) << "Invalid distribution pair: {" << entry.first << "," << entry.second << "}";
      ASSERT_TRUE(entry.second >= 0) << "Invalid distribution pair: {" << entry.first << "," << entry.second << "}";
      ASSERT_TRUE(entry.second < bulk.parallel_size()) << "Invalid distribution pair: {" << entry.first << "," << entry.second << "}";

      if(entry.second != bulk.parallel_rank()) {
        changeList.emplace_back(elem, entry.second);
      }
    }
  }

  bulk.change_entity_owner(changeList);
}

namespace simple_fields {

void ConstructedMesh::create_block_elements_and_nodes(stk::mesh::BulkData& bulk,
                                                      const ConstructedElementBlock& block,
                                                      const unsigned elemIdOffset)
{
  stk::mesh::Part* part = bulk.mesh_meta_data().get_part(block.name);
  STK_ThrowRequire(nullptr != part);

  size_t elementIndex = elemIdOffset;

  for(size_t elementCount = 0; elementCount < block.connectivityIndex.size(); ++elementIndex, ++elementCount) {
    stk::mesh::EntityIdVector elementNodeIds;
    for(size_t i=0; i<block.connectivityIndex[elementCount].size(); ++i) {
      int nodeIndex = block.connectivityIndex[elementCount][i];
      ASSERT_TRUE(nodeIndex > 0);
      stk::mesh::EntityId nodeId = m_nodeIds[nodeIndex - 1];
      elementNodeIds.push_back(nodeId);
    }
    stk::mesh::declare_element(bulk, *part, m_elemIds[elementIndex], elementNodeIds);
  }
}

void ConstructedMesh::verify_mesh_data(const stk::mesh::MetaData& meta)
{
  ASSERT_TRUE(meta.spatial_dimension() == 2 || meta.spatial_dimension() == 3);
  ASSERT_EQ(m_spatialDimension, meta.spatial_dimension());

  ASSERT_EQ(m_xCoords.size(), m_nodeIds.size());
  ASSERT_EQ(m_yCoords.size(), m_nodeIds.size());
  if(m_spatialDimension == 3) {
    ASSERT_EQ(m_zCoords.size(), m_nodeIds.size());
  }
}

void ConstructedMesh::populate_bulk_data(stk::mesh::BulkData& bulk)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  for(const ConstructedElementBlock& elemBlock : m_elemBlocks) {
    stk::mesh::Part& block = meta.declare_part_with_topology(elemBlock.name, elemBlock.topology);
    stk::io::put_io_part_attribute(block);
    meta.set_part_id(block, elemBlock.id);
  }

  stk::mesh::Field<double> & coordsField = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates", 1);
  stk::mesh::put_field_on_mesh(coordsField, meta.universal_part(), m_spatialDimension, nullptr);
  stk::io::set_field_output_type(coordsField, stk::io::FieldOutputType::VECTOR_3D);

  bulk.modification_begin();
  if(bulk.parallel_rank() == 0) {
    unsigned elemIdOffset = 0;
    for(unsigned i=0; i<m_elemBlocks.size(); ++i) {
      for(const std::vector<unsigned>& connectivity : m_elemBlocks[i].connectivityIndex) {
        ASSERT_EQ(m_elemBlocks[i].topology.num_nodes(), connectivity.size());
      }
      create_block_elements_and_nodes(bulk, m_elemBlocks[i], elemIdOffset);
      elemIdOffset += m_elemBlocks[i].connectivityIndex.size();
    }

    for(size_t nodeIndex=0; nodeIndex < m_nodeIds.size(); nodeIndex++)
    {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, m_nodeIds[nodeIndex]);
      ASSERT_TRUE(bulk.is_valid(node)) << "Invalid node id: " << m_nodeIds[nodeIndex];
      double * nodalCoords = stk::mesh::field_data(coordsField, node);

      nodalCoords[0] = m_xCoords[nodeIndex];
      nodalCoords[1] = m_yCoords[nodeIndex];

      if(m_spatialDimension == 3) {
        nodalCoords[2] = m_zCoords[nodeIndex];
      }
    }
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec changeList;

  if(bulk.parallel_rank() == 0 && bulk.parallel_size() > 1) {
    for(const stk::mesh::EntityIdProc& entry : m_elemDistribution) {
      stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, entry.first);
      ASSERT_TRUE(bulk.is_valid(elem)) << "Invalid distribution pair: {" << entry.first << "," << entry.second << "}";
      ASSERT_TRUE(entry.second >= 0) << "Invalid distribution pair: {" << entry.first << "," << entry.second << "}";
      ASSERT_TRUE(entry.second < bulk.parallel_size()) << "Invalid distribution pair: {" << entry.first << "," << entry.second << "}";

      if(entry.second != bulk.parallel_rank()) {
        changeList.emplace_back(elem, entry.second);
      }
    }
  }

  bulk.change_entity_owner(changeList);
}

} // namespace simple_fields

} // namespace unit_test_util
} // namespace stk

