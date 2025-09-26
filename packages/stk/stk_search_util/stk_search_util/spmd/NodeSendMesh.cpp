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
#include "stk_search_util/spmd/NodeSendMesh.hpp"
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"
#include "stk_topology/topology.hpp"                  // for operator<<, top...
#include "stk_mesh/base/Bucket.hpp"                   // for Bucket
#include <stk_mesh/base/BulkData.hpp>                 // for BulkData
#include <stk_mesh/base/Entity.hpp>                   // for Entity
#include <stk_mesh/base/FieldBase.hpp>                // for FieldBase
#include "stk_mesh/base/MetaData.hpp"                 // for MetaData
#include "stk_mesh/base/Selector.hpp"                 // for Selector
#include "stk_search/DistanceComparison.hpp"          // for distance_sq
#include "stk_search/SearchInterface.hpp"             // for FindParametricC...
#include "stk_util/parallel/ParallelReduce.hpp"       // for ReduceEnd, Reduce
#include "stk_util/util/ReportHandler.hpp"
#include <algorithm>                                  // for copy, sort, unique
#include <cmath>                                      // for sqrt
#include <memory>                                     // for make_shared
namespace stk::mesh { class Ghosting; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {
namespace spmd {

NodeSendMesh::NodeSendMesh(stk::mesh::BulkData* sendBulk, const stk::mesh::FieldBase* coordinateField,
                           const stk::mesh::PartVector& sendParts, const stk::ParallelMachine sendComm,
                           const double coarseSearchTolerance)
  : m_bulk(sendBulk)
  , m_meta(sendBulk != nullptr ? &sendBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_meshModified(false)
  , m_ghosting(nullptr)
  , m_syncCount(0)
  , m_coarseSearchTol(coarseSearchTolerance)
  , m_isInitialized(false)
  , m_extrapolateOption(stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG)
{
  STK_ThrowRequire(m_coordinateField->entity_rank() == stk::topology::NODE_RANK);

  if(sendBulk != nullptr) {
    m_findParametricCoords = std::make_shared<NodeParametricCoordsFinder>(*m_bulk, m_coordinateField);
    m_isInitialized = true;
  }
}

NodeSendMesh::NodeSendMesh(stk::mesh::BulkData* sendBulk, const stk::mesh::FieldBase* coordinateField,
                           const stk::mesh::PartVector& sendParts, const stk::mesh::Selector& activeSelector,
                           const stk::ParallelMachine sendComm, const double coarseSearchTolerance)
  : m_bulk(sendBulk)
  , m_meta(sendBulk != nullptr ? &sendBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(activeSelector)
  , m_meshModified(false)
  , m_ghosting(nullptr)
  , m_syncCount(0)
  , m_coarseSearchTol(coarseSearchTolerance)
  , m_isInitialized(false)
  , m_extrapolateOption(stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG)
{
  STK_ThrowRequire(m_coordinateField->entity_rank() == stk::topology::NODE_RANK);

  if(sendBulk != nullptr) {
    m_findParametricCoords = std::make_shared<NodeParametricCoordsFinder>(*m_bulk, m_coordinateField);
    m_isInitialized = true;
  }
}

void NodeSendMesh::consistency_check()
{
  STK_ThrowRequireMsg((m_meta != nullptr) && (m_bulk != nullptr),
                       "Meta and/or Bulk Data objects have not been set");

  for(const stk::mesh::Part* part : m_meshParts) {
    STK_ThrowRequireMsg((part != nullptr) && ((part->primary_entity_rank() <= stk::topology::ELEM_RANK) ||
                                              (part->primary_entity_rank() == stk::topology::INVALID_RANK)),
                        "All destination parts must be {ELEM_RANK} or {FACE_RANK} or {EDGE_RANK} or {NODE_RANK}");
  }
}

void NodeSendMesh::bounding_boxes(std::vector<BoundingBox>& v_box, bool includeGhosts) const
{
  STK_ThrowAssert(m_isInitialized);

  Point center;

  stk::mesh::EntityRank rank = stk::topology::NODE_RANK;
  stk::mesh::Selector selector = get_objects_selector(*m_meta, m_meshParts, &m_activeSelector, includeGhosts);

  stk::mesh::BucketVector const& buckets = m_bulk->get_buckets(rank, selector);

  const unsigned nDim = m_meta->spatial_dimension();

  for(auto&& ib : buckets) {
    stk::mesh::Bucket& b = *ib;

    for(auto node : b) {
      const double* coords = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, node));

      for(unsigned j = 0; j < nDim; ++j) {
        center[j] = coords[j];
      }

      EntityProc theIdent(make_entity_key_pair(m_bulk,node), m_bulk->parallel_rank());
      BoundingBox theBox(Sphere(center, m_coarseSearchTol), theIdent);
      v_box.push_back(theBox);
    }
  }
  std::sort(v_box.begin(), v_box.end(),
            [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
}

void NodeSendMesh::find_parametric_coords(const EntityKey& k, const std::vector<double>& toCoords,
                                          std::vector<double>& parametricCoords,
                                          double& parametricDistance,
                                          bool& isWithinParametricTolerance) const
{
  STK_ThrowAssert(m_isInitialized);
  m_findParametricCoords->find_parametric_coords(k, toCoords,
                                                 parametricCoords,
                                                 parametricDistance,
                                                 isWithinParametricTolerance);
}

bool NodeSendMesh::modify_search_outside_parametric_tolerance(const EntityKey& /*k*/,
                                                              const std::vector<double>& /*toCoords*/,
                                                              std::vector<double>& /*parametricCoords*/,
                                                              double& /*geometricDistanceSquared*/,
                                                              bool& /*isWithinGeometricTolerance*/) const
{
  return false;
}

double
NodeSendMesh::get_closest_geometric_distance_squared(const EntityKey& k, const std::vector<double>& toCoords) const
{
  unsigned numDimension = m_meta->spatial_dimension();
  const stk::mesh::Entity node = k;
  const double* entityCoords = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, node));
  return stk::search::distance_sq(numDimension, entityCoords, toCoords.data());
}

void NodeSendMesh::initialize()
{
  m_isInitialized = m_bulk != nullptr;

  if(m_bulk != nullptr && !m_findParametricCoords) {
    m_findParametricCoords = std::make_shared<NodeParametricCoordsFinder>(*m_bulk, m_coordinateField);
    m_isInitialized = true;
  }

  consistency_check();
}

void NodeSendMesh::update_ghosting(const EntityProcVec& entity_keys, const std::string& suffix)
{
  STK_ThrowAssert(m_isInitialized);
  m_ghostingMap.clear();

  for(auto key_proc : entity_keys) {
    // convert from EntityProc based on EntityKey to EntityProc based on raw Entity.
    const unsigned proc = key_proc.proc();
    const stk::mesh::Entity e = key_proc.id();
    const stk::mesh::EntityProc ep(e, proc);

    m_ghostingMap.push_back(ep);
  }

  unsigned s = !m_ghostingMap.empty();

  stk::all_reduce(m_comm, stk::ReduceSum<1>(&s));

  if(s) {
    std::sort(m_ghostingMap.begin(), m_ghostingMap.end());
    m_ghostingMap.erase(unique(m_ghostingMap.begin(), m_ghostingMap.end()), m_ghostingMap.end());

    std::string theGhostName = "STK_transfer_ghosting" + suffix;

    if(m_ghosting == nullptr) {
      const bool iStartedMod = m_bulk->modification_begin();

      m_ghosting = &m_bulk->create_ghosting(theGhostName);

      if (iStartedMod) m_bulk->modification_end();
    }

    const bool iStartedMod = m_bulk->modification_begin();
    m_bulk->change_ghosting(*m_ghosting, m_ghostingMap);
    if (iStartedMod) m_bulk->modification_end();

    m_meshModified = true;
  }
}

void NodeSendMesh::update_ghosted_key(EntityKey& k)
{
  stk::mesh::EntityKey key = k;
  stk::mesh::Entity entity = k;

  stk::mesh::Entity keyEntity = m_bulk->get_entity(key);

  if(keyEntity != entity) {
    k.m_entity = keyEntity;
  }
}

void NodeSendMesh::destroy_ghosting()
{
  STK_ThrowAssert(m_isInitialized);
  if(m_ghosting != nullptr) {
    m_bulk->destroy_ghosting(*m_ghosting);
  }
}

bool NodeSendMesh::is_valid(const EntityKey& k) const {
  STK_ThrowAssert(m_isInitialized);
  return m_bulk->is_valid(k);
}

double NodeSendMesh::get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& toCoords) const
{
  return get_distance_from_centroid(k, toCoords);
}

double NodeSendMesh::get_distance_from_centroid(const EntityKey& k, const std::vector<double>& toCoords) const
{
  double distanceSquared = get_distance_squared_from_centroid(k, toCoords);
  return std::sqrt(distanceSquared);
}

double NodeSendMesh::get_distance_squared_from_centroid(const EntityKey& k, const std::vector<double>& toCoords) const
{
  STK_ThrowAssert(m_isInitialized);
  const stk::mesh::Entity e = k;
  const unsigned nDim = m_meta->spatial_dimension();

  const double* coor = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, e));
  return stk::search::distance_sq(nDim, coor, toCoords.data());
}

stk::mesh::EntityId NodeSendMesh::id(const EntityKey& k) const
{
  STK_ThrowAssert(m_isInitialized);
  return m_bulk->identifier(k);
}

void NodeSendMesh::centroid(const EntityKey& k, std::vector<double>& centroid) const
{
  const stk::mesh::Entity node = k;
  const double* centroidCoord = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, node));
  const unsigned nDim = m_meta->spatial_dimension();
  centroid.assign(centroidCoord, centroidCoord + nDim);
}

void NodeSendMesh::coordinates(const EntityKey& k, std::vector<double>& coords) const
{
  STK_ThrowAssert(m_isInitialized);
  const stk::mesh::Entity node = k;
  const double* data = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, node));
  const unsigned nDim = m_meta->spatial_dimension();
  coords.assign(data, data+nDim);
}

void NodeSendMesh::post_mesh_modification_event()
{
  STK_ThrowAssert(m_isInitialized);
  m_syncCount = m_bulk->synchronized_count();
}

std::vector<std::string> NodeSendMesh::get_part_membership(const EntityKey& k) const
{
  const stk::mesh::Entity elem = k;
  return stk::search::get_part_membership(*m_bulk, elem, m_meshParts);
}

} // namespace spmd
} // namespace search
} // namespace stk
