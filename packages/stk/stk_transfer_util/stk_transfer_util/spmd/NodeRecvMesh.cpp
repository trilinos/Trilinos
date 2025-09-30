// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_transfer_util/spmd/NodeRecvMesh.hpp"
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"
#include <stk_mesh/base/Entity.hpp>         // for Entity
#include <stk_mesh/base/FieldBase.hpp>      // for FieldBase
#include "stk_mesh/base/MetaData.hpp"       // for MetaData
#include "stk_mesh/base/Selector.hpp"       // for Selector
#include "stk_mesh/base/Bucket.hpp"         // for Bucket
#include <stk_mesh/base/BulkData.hpp>       // for BulkData
#include "stk_topology/topology.hpp"        // for topology, operator<<, top...
#include "stk_util/util/ReportHandler.hpp"  // for eval_test_condition, STK_...

#include <algorithm>                        // for max, sort, fill
#include <cstddef>                          // for size_t
#include <memory>                           // for shared_ptr, __shared_ptr_...
#include <type_traits>

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {
namespace spmd {

NodeRecvMesh::NodeRecvMesh(stk::mesh::BulkData* recvBulk,
                           const stk::mesh::FieldBase* coordinateField,
                           const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                           const stk::mesh::PartVector& recvParts,
                           const stk::ParallelMachine recvComm,
                           double parametricTolerance, double geometricTolerance)
  : m_bulk(recvBulk)
  , m_meta(recvBulk != nullptr ? &recvBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_fieldSpecs(fieldSpecs)
  , m_fieldVec(stk::transfer::get_fields(m_meta, m_fieldSpecs, stk::topology::NODE_RANK))
  , m_meshParts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
  , m_searchMesh(m_bulk, m_coordinateField, m_meshParts, m_activeSelector,
                 m_comm, m_parametricTolerance, m_searchTolerance)
{
  initialize();
}

NodeRecvMesh::NodeRecvMesh(stk::mesh::BulkData* recvBulk,
                           const stk::mesh::FieldBase* coordinateField,
                           const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                           const stk::mesh::PartVector& recvParts,
                           const stk::mesh::Selector& activeSelector,
                           const stk::ParallelMachine recvComm,
                           double parametricTolerance, double geometricTolerance)
  : m_bulk(recvBulk)
  , m_meta(recvBulk != nullptr ? &recvBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_fieldSpecs(fieldSpecs)
  , m_fieldVec(stk::transfer::get_fields(m_meta, m_fieldSpecs, stk::topology::NODE_RANK))
  , m_meshParts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(activeSelector)
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
  , m_searchMesh(m_bulk, m_coordinateField, m_meshParts, m_activeSelector,
                 m_comm, m_parametricTolerance, m_searchTolerance)
{
  initialize();
}

void NodeRecvMesh::consistency_check()
{
  static_assert( std::is_same_v<TransferBaseClass::Entity,        SearchBaseClass::Entity>         == true );
  static_assert( std::is_same_v<TransferBaseClass::EntityVec,     SearchBaseClass::EntityVec>      == true );
  static_assert( std::is_same_v<TransferBaseClass::EntityKey,     SearchBaseClass::EntityKey>      == true );
  static_assert( std::is_same_v<TransferBaseClass::EntityProc,    SearchBaseClass::EntityProc>     == true );
  static_assert( std::is_same_v<TransferBaseClass::EntityProcVec, SearchBaseClass::EntityProcVec>  == true );
  static_assert( std::is_same_v<TransferBaseClass::Point,         SearchBaseClass::Point>          == true );
  static_assert( std::is_same_v<TransferBaseClass::Box,           SearchBaseClass::Box>            == true );
  static_assert( std::is_same_v<TransferBaseClass::Sphere,        SearchBaseClass::Sphere>         == true );
  static_assert( std::is_same_v<TransferBaseClass::BoundingBox,   SearchBaseClass::BoundingBox>    == true );

  for(const stk::mesh::Part* part : m_meshParts) {
    STK_ThrowRequireMsg((part != nullptr) && ((part->primary_entity_rank() <= stk::topology::ELEM_RANK) ||
                                              (part->primary_entity_rank() == stk::topology::INVALID_RANK)),
                        "All destination parts must be {ELEM_RANK} or {FACE_RANK} or {EDGE_RANK} or {NODE_RANK}");
  }
}

void NodeRecvMesh::initialize()
{
  m_isInitialized = (m_bulk != nullptr);
  m_searchMesh.initialize();

  consistency_check();
}

void NodeRecvMesh::initialize(stk::mesh::BulkData* recvBulk)
{
  if(nullptr != recvBulk) {
    m_bulk = recvBulk;
    m_meta = &recvBulk->mesh_meta_data();
    m_fieldVec = stk::transfer::get_fields(m_meta, m_fieldSpecs, stk::topology::NODE_RANK);
    m_isInitialized = true;

    m_searchMesh.initialize(m_bulk);

    consistency_check();
  }
}

unsigned NodeRecvMesh::spatial_dimension() const
{
  STK_ThrowAssert(m_isInitialized);
  return m_meta->spatial_dimension();
}

void NodeRecvMesh::bounding_boxes(std::vector<BoundingBox>& v_box) const
{
  STK_ThrowAssert(m_isInitialized);
  m_searchMesh.bounding_boxes(v_box);
}

void NodeRecvMesh::coordinates(const EntityKey& k, std::vector<double>& coords) const
{
  STK_ThrowAssert(m_isInitialized);
  m_searchMesh.coordinates(k, coords);
}

stk::mesh::EntityId NodeRecvMesh::id(const EntityKey& k) const
{
  STK_ThrowAssert(m_isInitialized);
  return m_searchMesh.id(k);
}

double NodeRecvMesh::get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& toCoords) const
{
  STK_ThrowAssert(m_isInitialized);
  return m_searchMesh.get_distance_from_nearest_node(k, toCoords);
}

void NodeRecvMesh::centroid(const EntityKey& k, std::vector<double>& centroid) const
{
  STK_ThrowAssert(m_isInitialized);
  m_searchMesh.centroid(k, centroid);
}

void NodeRecvMesh::fill_entity_keys(const stk::mesh::EntityKeyVector& rangeEntities,
                                            std::vector<EntityKey>& elementEntityKeys)
{
  STK_ThrowAssert(m_isInitialized);
  m_searchMesh.fill_entity_keys(rangeEntities, elementEntityKeys);
}

double* NodeRecvMesh::value(const EntityKey& k, const unsigned i) const
{
  STK_ThrowAssert(m_isInitialized);

  const stk::mesh::FieldBase* field = m_fieldVec[i].field;

  field->sync_to_host();
  field->modify_on_host();

  double* data = static_cast<double *>(stk::mesh::field_data(*field, k));
  return data;
}

unsigned NodeRecvMesh::value_size(const EntityKey& k, const unsigned i) const
{
  STK_ThrowAssert(m_isInitialized);

  return stk::mesh::field_scalars_per_entity(*m_fieldVec[i].field, k);
}

unsigned NodeRecvMesh::num_values(const EntityKey& /*e*/) const { return m_fieldVec.size(); }

unsigned NodeRecvMesh::max_num_values() const { return m_fieldVec.size(); }

unsigned NodeRecvMesh::value_key(const EntityKey& /*k*/, const unsigned i) const
{
  return i;
}

void NodeRecvMesh::update_values()
{
  STK_ThrowAssert(m_isInitialized);

  for(auto field : m_fieldVec) {
    field.field->sync_to_host();
    field.field->modify_on_host();
  }

  std::vector<const stk::mesh::FieldBase*> fields = stk::transfer::extract_field_pointers(m_fieldVec);
  stk::mesh::communicate_field_data(*m_bulk, fields);
}

unsigned NodeRecvMesh::get_index(const unsigned i) const
{
  return m_fieldVec[i].index;
}

} // namespace spmd
} // namespace transfer
} // namespace stk

