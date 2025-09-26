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
#include "stk_search_util/spmd/ElementSendMesh.hpp"
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"
#include "stk_search_util/PointProjection.hpp"
#include "stk_mesh/base/Selector.hpp"                 // for Selector
#include "stk_search/SearchInterface.hpp"             // for ExternalPointHa...
#include "stk_util/parallel/ParallelReduce.hpp"       // for ReduceEnd, Reduce
#include <stk_mesh/base/Entity.hpp>                   // for Entity
#include "stk_util/util/ReportHandler.hpp"
#include "stk_mesh/base/Bucket.hpp"                   // for Bucket
#include <stk_mesh/base/BulkData.hpp>                 // for BulkData, commu...
#include <stk_mesh/base/FieldBase.hpp>                // for FieldBase
#include "stk_mesh/base/MetaData.hpp"                 // for MetaData
#include "stk_search/DistanceComparison.hpp"          // for distance_sq
#include "stk_topology/topology.hpp"                  // for operator<<, top...
#include <algorithm>                                  // for max, sort, min
#include <cmath>                                      // for sqrt
#include <memory>                                     // for shared_ptr, __s...
#include <limits>
namespace stk::mesh { class Ghosting; }
namespace stk::mesh { class Part; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {
namespace spmd {

namespace impl {
void fill_bounding_box(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                       const stk::mesh::FieldBase* coordField,
                       stk::search::Point<double>& minCorner,
                       stk::search::Point<double>& maxCorner)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const unsigned nDim = meta.spatial_dimension();

  STK_ThrowRequireMsg(bulk.is_valid(entity), "Invalid entity: " << bulk.entity_key(entity));

  double doubleMax = std::numeric_limits<double>::max();

  for(unsigned j = 0; j < nDim; ++j) {
    minCorner[j] = +doubleMax;
    maxCorner[j] = -doubleMax;
  }

  stk::mesh::Entity const* nodes = bulk.begin_nodes(entity);
  unsigned numNodes = bulk.num_nodes(entity);
  for(unsigned ni = 0; ni < numNodes; ++ni) {
    stk::mesh::Entity node = nodes[ni];

    const double* coords = static_cast<const double *>(stk::mesh::field_data(*coordField, node));

    for(unsigned j = 0; j < nDim; ++j) {
      minCorner[j] = std::min(minCorner[j], coords[j]);
      maxCorner[j] = std::max(maxCorner[j], coords[j]);
    }
  }
}
}

ElementSendMesh::ElementSendMesh(stk::mesh::BulkData* sendBulk, const stk::mesh::FieldBase* coordinateField,
                                 const stk::mesh::EntityRank sendEntityRank, const stk::mesh::PartVector& sendParts,
                                 const stk::ParallelMachine sendComm,
                                 std::shared_ptr<FindParametricCoordsInterface> findParametricCoords,
                                 std::shared_ptr<HandleExternalPointInterface> externalPointHandler,
                                 std::shared_ptr<MasterElementProviderInterface> masterElemProvider)
  : m_bulk(sendBulk)
  , m_meta(sendBulk != nullptr ? &sendBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_searchEntityRank(sendEntityRank)
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_meshModified(false)
  , m_ghosting(nullptr)
  , m_findParametricCoords(findParametricCoords)
  , m_externalPointHandler(externalPointHandler)
  , m_masterElementProvider(masterElemProvider)
  , m_syncCount(0)
  , m_isInitialized(false)
  , m_extrapolateOption(stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG)
{
}

ElementSendMesh::ElementSendMesh(stk::mesh::BulkData* sendBulk, const stk::mesh::FieldBase* coordinateField,
                                 const stk::mesh::EntityRank sendEntityRank, const stk::mesh::PartVector& sendParts,
                                 const stk::mesh::Selector& activeSelector,
                                 const stk::ParallelMachine sendComm,
                                 std::shared_ptr<FindParametricCoordsInterface> findParametricCoords,
                                 std::shared_ptr<HandleExternalPointInterface> externalPointHandler,
                                 std::shared_ptr<MasterElementProviderInterface> masterElemProvider)
  : m_bulk(sendBulk)
  , m_meta(sendBulk != nullptr ? &sendBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_searchEntityRank(sendEntityRank)
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(activeSelector)
  , m_meshModified(false)
  , m_ghosting(nullptr)
  , m_findParametricCoords(findParametricCoords)
  , m_externalPointHandler(externalPointHandler)
  , m_masterElementProvider(masterElemProvider)
  , m_syncCount(0)
  , m_isInitialized(false)
  , m_extrapolateOption(stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG)
{
}

ElementSendMesh::ElementSendMesh(stk::mesh::BulkData* sendBulk, const stk::mesh::FieldBase* coordinateField,
                                 const stk::mesh::EntityRank sendEntityRank, const stk::mesh::PartVector& sendParts,
                                 const stk::ParallelMachine sendComm,
                                 std::shared_ptr<FindParametricCoordsInterface> findParametricCoords,
                                 std::shared_ptr<HandleExternalPointInterface> externalPointHandler)
  : m_bulk(sendBulk)
  , m_meta(sendBulk != nullptr ? &sendBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_searchEntityRank(sendEntityRank)
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_meshModified(false)
  , m_ghosting(nullptr)
  , m_findParametricCoords(findParametricCoords)
  , m_externalPointHandler(externalPointHandler)
  , m_syncCount(0)
  , m_isInitialized(false)
  , m_extrapolateOption(stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG)
{
}

ElementSendMesh::ElementSendMesh(stk::mesh::BulkData* sendBulk, const stk::mesh::FieldBase* coordinateField,
                                 const stk::mesh::EntityRank sendEntityRank, const stk::mesh::PartVector& sendParts,
                                 const stk::mesh::Selector& activeSelector,
                                 const stk::ParallelMachine sendComm,
                                 std::shared_ptr<FindParametricCoordsInterface> findParametricCoords,
                                 std::shared_ptr<HandleExternalPointInterface> externalPointHandler)
  : m_bulk(sendBulk)
  , m_meta(sendBulk != nullptr ? &sendBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_searchEntityRank(sendEntityRank)
  , m_meshParts(sendParts)
  , m_comm(sendComm)
  , m_activeSelector(activeSelector)
  , m_meshModified(false)
  , m_ghosting(nullptr)
  , m_findParametricCoords(findParametricCoords)
  , m_externalPointHandler(externalPointHandler)
  , m_syncCount(0)
  , m_isInitialized(false)
  , m_extrapolateOption(stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG)
{
}

void ElementSendMesh::initialize()
{
  m_isInitialized = m_bulk != nullptr;

  consistency_check();
}

void ElementSendMesh::consistency_check()
{
  STK_ThrowRequireMsg((m_meta != nullptr) && (m_bulk != nullptr),
                       "Meta and/or Bulk Data objects have not been set");

  STK_ThrowRequireMsg(
      (m_searchEntityRank == stk::topology::ELEM_RANK) || (m_searchEntityRank == stk::topology::FACE_RANK) ||
      (m_searchEntityRank == stk::topology::EDGE_RANK && m_meta->spatial_dimension() == 2),
      "Input object type must be {ELEM_RANK} or {FACE_RANK} or {EDGE_RANK with dimension 2}: " << m_searchEntityRank);

  for(const stk::mesh::Part* part : m_meshParts) {
    STK_ThrowRequireMsg(part_has_proper_entity_rank(part),
                    "All source parts must be {ELEM_RANK} or {FACE_RANK} or {EDGE_RANK with dimension 2}");
  }
}

void ElementSendMesh::bounding_boxes(std::vector<BoundingBox>& v_box, bool includeGhosts) const
{
  STK_ThrowAssert(m_isInitialized);

  Point minCorner, maxCorner;

  stk::mesh::EntityRank rank = get_objects_rank(m_searchEntityRank, m_meshParts);
  stk::mesh::Selector selector = get_objects_selector(*m_meta, m_meshParts, &m_activeSelector, includeGhosts);

  stk::mesh::BucketVector const& buckets = m_bulk->get_buckets(rank, selector);

  for(auto&& ib : buckets) {
    stk::mesh::Bucket& b = *ib;

    for(auto entity : b) {
      impl::fill_bounding_box(*m_bulk, entity, m_coordinateField, minCorner, maxCorner);

      EntityKey key(entity, m_bulk->entity_key(entity));
      EntityProc theIdent(key, m_bulk->parallel_rank());
      BoundingBox theBox(Box(minCorner, maxCorner), theIdent);
      v_box.push_back(theBox);
    }
  }
  std::sort(v_box.begin(), v_box.end(),
            [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
}


void ElementSendMesh::find_parametric_coords(const EntityKey& k, const std::vector<double>& toCoords,
                                             std::vector<double>& parametricCoords,
                                             double& parametricDistance,
                                             bool& isWithinParametricTolerance) const
{
  m_findParametricCoords->find_parametric_coords(k, toCoords,
                                                 parametricCoords,
                                                 parametricDistance,
                                                 isWithinParametricTolerance);
}

bool ElementSendMesh::modify_search_outside_parametric_tolerance(const EntityKey& k,
                                                                 const std::vector<double>& toCoords,
                                                                 std::vector<double>& parametricCoords,
                                                                 double& geometricDistanceSquared,
                                                                 bool& isWithinGeometricTolerance) const
{
  return m_externalPointHandler->handle_point(k, toCoords, parametricCoords, geometricDistanceSquared, isWithinGeometricTolerance);
}

double ElementSendMesh::get_closest_geometric_distance_squared(const EntityKey& k, const std::vector<double>& toCoords) const
{
  STK_ThrowAssert(m_isInitialized);
  stk::mesh::Entity e = k;

  if(m_masterElementProvider) {
    stk::search::ProjectionData data(*m_bulk, m_masterElementProvider, toCoords, *m_coordinateField);
    stk::search::ProjectionResult projectionResult;
    if(k.rank() == stk::topology::ELEM_RANK) {
      stk::search::project_to_closest_side(data, e, projectionResult);
    }
    else {
      stk::search::project_to_entity(data, e, projectionResult);
    }
    STK_ThrowRequire(projectionResult.doneProjection);
    return projectionResult.geometricDistanceSquared;
  }

  return distance_squared_from_nearest_entity_node(*m_bulk, e, m_coordinateField, toCoords);
}

void ElementSendMesh::update_ghosting(const EntityProcVec& entity_keys, const std::string& suffix)
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

    std::vector<const stk::mesh::FieldBase*> fields;
    fields.push_back(m_coordinateField);
    stk::mesh::communicate_field_data(*m_bulk, fields);

    m_meshModified = true;
  }
}

void ElementSendMesh::update_ghosted_key(EntityKey& k)
{
  stk::mesh::EntityKey key = k;
  stk::mesh::Entity entity = k;

  stk::mesh::Entity keyEntity = m_bulk->get_entity(key);

  if(keyEntity != entity) {
    k.m_entity = keyEntity;
  }
}

void ElementSendMesh::destroy_ghosting()
{
  STK_ThrowAssert(m_isInitialized);
  if(m_ghosting != nullptr) {
    m_bulk->destroy_ghosting(*m_ghosting);
  }
}

bool ElementSendMesh::is_valid(const EntityKey& k) const {
  STK_ThrowAssert(m_isInitialized);
  return m_bulk->is_valid(k);
}

double ElementSendMesh::get_distance_from_centroid(const EntityKey& k, const std::vector<double>& toCoords) const
{
  double distanceSquared = get_distance_squared_from_centroid(k, toCoords);
  return std::sqrt(distanceSquared);
}

double ElementSendMesh::get_distance_squared_from_centroid(const EntityKey& k, const std::vector<double>& toCoords) const
{
  STK_ThrowAssert(m_isInitialized);
  const stk::mesh::Entity e = k;
  double distanceSquared = std::numeric_limits<double>::max();
  const unsigned nDim = m_meta->spatial_dimension();

  if(m_coordinateField->entity_rank() == stk::topology::NODE_RANK) {
    double coordVector[3] = {0.0};
    determine_centroid(nDim, e, *m_coordinateField, coordVector);
    distanceSquared = stk::search::distance_sq(nDim, coordVector, toCoords.data());
  }
  else if(m_coordinateField->entity_rank() == stk::topology::ELEM_RANK) {
    const double* coor = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, e));
    distanceSquared = stk::search::distance_sq(nDim, coor, toCoords.data());
  }

  return distanceSquared;
}

double ElementSendMesh::get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& toCoords) const
{
  STK_ThrowAssert(m_isInitialized);
  const stk::mesh::Entity e = k;
  return distance_from_nearest_entity_node(*m_bulk, e, m_coordinateField, toCoords);
}

stk::mesh::EntityId ElementSendMesh::id(const EntityKey& k) const
{
  STK_ThrowAssert(m_isInitialized);
  stk::mesh::Entity e = k;
  return m_bulk->identifier(e);
}

void ElementSendMesh::centroid(const EntityKey& k, std::vector<double>& centroid) const
{
  coordinates(k, centroid);
}

void ElementSendMesh::coordinates(const EntityKey& k, std::vector<double>& coords) const
{
  STK_ThrowAssert(m_isInitialized);
  const stk::mesh::Entity elem = k;
  const unsigned nDim = m_meta->spatial_dimension();
  if(m_coordinateField->entity_rank() == stk::topology::NODE_RANK) {
    determine_centroid(nDim, elem, *m_coordinateField, coords);
  }
  else if(m_coordinateField->entity_rank() == stk::topology::ELEM_RANK) {
    const double* coor = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, elem));
    coords.assign(coor, coor+nDim);
  }
}

void ElementSendMesh::post_mesh_modification_event()
{
  m_syncCount = m_bulk->synchronized_count();
}

std::vector<std::string> ElementSendMesh::get_part_membership(const EntityKey& k) const
{
  const stk::mesh::Entity elem = k;
  return stk::search::get_part_membership(*m_bulk, elem, m_meshParts);
}

} // namespace spmd
} // namespace search
} // namespace stk

