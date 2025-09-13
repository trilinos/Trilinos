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
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"
#include "stk_mesh/base/Bucket.hpp"         // for Bucket
#include <stk_mesh/base/BulkData.hpp>       // for BulkData
#include <stk_mesh/base/Entity.hpp>         // for Entity
#include <stk_mesh/base/FieldBase.hpp>      // for FieldBase
#include "stk_mesh/base/MetaData.hpp"       // for MetaData
#include "stk_mesh/base/Part.hpp"                         // for Part
#include "stk_mesh/base/Selector.hpp"       // for Selector
#include "stk_topology/topology.hpp"        // for topology, operator<<, top...
#include "stk_transfer_util/PointInterpolation.hpp"
#include "stk_transfer_util/spmd/NodeSendMesh.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for eval_test_condition, STK_...

#include <algorithm>                                      // for max
#include <cstddef>                                        // for size_t
#include <memory>                                         // for shared_ptr
namespace stk::mesh { class Ghosting; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {
namespace spmd {

NodeSendMesh::NodeSendMesh(stk::mesh::BulkData* sendBulk,
                           const stk::mesh::FieldBase* coordinateField,
                           const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                           const stk::mesh::PartVector& sendParts,
                           const stk::ParallelMachine sendComm,
                           const double coarseSearchTolerance)
  : m_bulk(sendBulk)
  , m_meta(m_bulk != nullptr ? &m_bulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_fieldSpecs(fieldSpecs)
  , m_fieldVec(stk::transfer::get_fields(m_meta, m_fieldSpecs, stk::topology::NODE_RANK))
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_searchMesh(m_bulk, m_coordinateField, m_meshParts, m_activeSelector, m_comm, coarseSearchTolerance)
{
  initialize();
}

NodeSendMesh::NodeSendMesh(stk::mesh::BulkData* sendBulk,
                           const stk::mesh::FieldBase* coordinateField,
                           const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                           const stk::mesh::PartVector& sendParts,
                           const stk::ParallelMachine sendComm,
                           std::shared_ptr<stk::transfer::InterpolateFieldsInterface> interpolateFields,
                           const double coarseSearchTolerance)
  : m_bulk(sendBulk)
  , m_meta(m_bulk != nullptr ? &m_bulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_fieldSpecs(fieldSpecs)
  , m_fieldVec(stk::transfer::get_fields(m_meta, m_fieldSpecs, stk::topology::NODE_RANK))
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_interpolateFields(interpolateFields)
  , m_searchMesh(m_bulk, m_coordinateField, m_meshParts, m_activeSelector, m_comm, coarseSearchTolerance)
{
  initialize();
}

NodeSendMesh::NodeSendMesh(stk::mesh::BulkData* sendBulk,
                           const stk::mesh::FieldBase* coordinateField,
                           const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                           const stk::mesh::PartVector& sendParts,
                           const stk::mesh::Selector& activeSelector,
                           const stk::ParallelMachine sendComm,
                           const double coarseSearchTolerance)
  : m_bulk(sendBulk)
  , m_meta(m_bulk != nullptr ? &m_bulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_fieldSpecs(fieldSpecs)
  , m_fieldVec(stk::transfer::get_fields(m_meta, m_fieldSpecs, stk::topology::NODE_RANK))
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(activeSelector)
  , m_searchMesh(m_bulk, m_coordinateField, m_meshParts, m_activeSelector, m_comm, coarseSearchTolerance)
{
  initialize();
}

NodeSendMesh::NodeSendMesh(stk::mesh::BulkData* sendBulk,
                           const stk::mesh::FieldBase* coordinateField,
                           const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                           const stk::mesh::PartVector& sendParts,
                           const stk::mesh::Selector& activeSelector,
                           const stk::ParallelMachine sendComm,
                           std::shared_ptr<stk::transfer::InterpolateFieldsInterface> interpolateFields,
                           const double coarseSearchTolerance)
  : m_bulk(sendBulk)
  , m_meta(m_bulk != nullptr ? &m_bulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_fieldSpecs(fieldSpecs)
  , m_fieldVec(stk::transfer::get_fields(m_meta, m_fieldSpecs, stk::topology::NODE_RANK))
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(activeSelector)
  , m_interpolateFields(interpolateFields)
  , m_searchMesh(m_bulk, m_coordinateField, m_meshParts, m_activeSelector, m_comm, coarseSearchTolerance)
{
  initialize();
}

void NodeSendMesh::initialize()
{
  STK_ThrowRequire(m_coordinateField->entity_rank() == stk::topology::NODE_RANK);

  m_isInitialized = (m_bulk != nullptr);
  m_searchMesh.initialize();

  if(m_bulk != nullptr) {
    if(!m_interpolateFields) {
      m_interpolateFields = std::make_shared<CopyFieldInterpolator>(*m_bulk);
    }
  }

  consistency_check();
}

void NodeSendMesh::consistency_check()
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

void NodeSendMesh::bounding_boxes(std::vector<TransferBaseClass::BoundingBox>& v_box, bool includeGhosts) const
{
  STK_ThrowAssert(m_isInitialized);
  m_searchMesh.bounding_boxes(v_box, includeGhosts);
}

std::vector<std::string> NodeSendMesh::get_part_membership(const EntityKey& k) const
{
  return m_searchMesh.get_part_membership(k);
}

std::string NodeSendMesh::get_inspector_delimiter() const
{
  return m_searchMesh.get_inspector_delimiter();
}

void NodeSendMesh::find_parametric_coords(const EntityKey& k,
                                          const std::vector<double>& evalPoint,
                                          std::vector<double>& parametricCoords,
                                          double& parametricDistance,
                                          bool& isWithinParametricTolerance) const
{
  m_searchMesh.find_parametric_coords(k, evalPoint, parametricCoords, parametricDistance, isWithinParametricTolerance);
}

bool NodeSendMesh::modify_search_outside_parametric_tolerance(const EntityKey& k,
                                                              const std::vector<double>& evalPoint,
                                                              std::vector<double>& parametricCoords,
                                                              double& geometricDistanceSquared,
                                                              bool& isWithinGeometricTolerance) const
{
  return m_searchMesh.modify_search_outside_parametric_tolerance(k, evalPoint, parametricCoords,
                                                                 geometricDistanceSquared, isWithinGeometricTolerance);
}

double NodeSendMesh::get_closest_geometric_distance_squared(const EntityKey& k, const std::vector<double>& evalPoint) const
{
  return m_searchMesh.get_closest_geometric_distance_squared(k, evalPoint);
}

void NodeSendMesh::interpolate_fields(const EntityKey& k, std::vector<double>& evalPoint,
                                         std::vector<double>& sendParametricCoords, InterpolationData& data) const
{
  m_interpolateFields->interpolate_fields(k, evalPoint, sendParametricCoords, data);
}

void NodeSendMesh::update_ghosting(const EntityProcVec& entity_keys, const std::string& suffix)
{
  std::string searchSuffix = suffix;
  for(unsigned i = 0; i != m_fieldVec.size(); ++i) {
    searchSuffix += "_" + m_fieldVec[i].field->name();
  }

  m_searchMesh.update_ghosting(entity_keys, searchSuffix);
}

void NodeSendMesh::update_ghosted_key(EntityKey& k)
{
  m_searchMesh.update_ghosted_key(k);
}

void NodeSendMesh::destroy_ghosting()
{
  m_searchMesh.destroy_ghosting();
}

void NodeSendMesh::update_values()
{
  STK_ThrowAssert(m_isInitialized);

  for(auto field : m_fieldVec) {
    field.field->sync_to_host();
    field.field->modify_on_host();
  }

  stk::mesh::Ghosting* ghosting = m_searchMesh.get_ghosting();
  if(nullptr != ghosting) {
    std::vector<const stk::mesh::FieldBase*> fields = stk::transfer::extract_field_pointers(m_fieldVec);
    if(m_searchMesh.is_mesh_modified()) {
      // Copy coordinates to the newly ghosted nodes
      m_searchMesh.set_mesh_modified(false);
      fields.push_back(m_coordinateField);
    }
    stk::mesh::communicate_field_data(*m_bulk, fields);
  }
}

bool NodeSendMesh::is_valid(const EntityKey& k) const { return m_bulk->is_valid(k); }

double NodeSendMesh::get_distance_from_centroid(const EntityKey& k, const std::vector<double>& evalPoint) const
{
  return m_searchMesh.get_distance_from_centroid(k, evalPoint);
}

double NodeSendMesh::get_distance_squared_from_centroid(const EntityKey& k, const std::vector<double>& evalPoint) const
{
  return m_searchMesh.get_distance_squared_from_centroid(k, evalPoint);
}

double NodeSendMesh::get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& evalPoint) const
{
  return m_searchMesh.get_distance_from_nearest_node(k, evalPoint);
}

stk::mesh::EntityId NodeSendMesh::id(const EntityKey& k) const
{
  return m_bulk->identifier(k);
}

void NodeSendMesh::centroid(const EntityKey& k, std::vector<double>& centroid) const
{
  m_searchMesh.centroid(k, centroid);
}

void NodeSendMesh::coordinates(const EntityKey& k, std::vector<double>& coords) const
{
  m_searchMesh.coordinates(k, coords);
}

void NodeSendMesh::post_mesh_modification_event()
{
  m_searchMesh.post_mesh_modification_event();
}

void NodeSendMesh::set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy option)
{
  m_searchMesh.set_extrapolate_option(option);
}

stk::search::ObjectOutsideDomainPolicy NodeSendMesh::get_extrapolate_option() const
{
  return m_searchMesh.get_extrapolate_option();
}


} // namespace spmd
} // namespace transfer
} // namespace stk

