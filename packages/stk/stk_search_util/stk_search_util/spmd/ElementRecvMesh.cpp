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
#include "stk_search_util/spmd/ElementRecvMesh.hpp"
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

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {
namespace spmd {

ElementRecvMesh::ElementRecvMesh(stk::mesh::BulkData* recvBulk,
                                 const stk::mesh::FieldBase* coordinateField,
                                 const stk::mesh::EntityRank searchEntityRank,
                                 const stk::mesh::PartVector& recvParts,
                                 const stk::ParallelMachine recvComm,
                                 std::shared_ptr<PointEvaluatorInterface> pointEvaluator,
                                 double parametricTolerance, double geometricTolerance)
  : m_bulk(recvBulk)
  , m_meta(recvBulk != nullptr ? &recvBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_searchEntityRank(searchEntityRank)
  , m_meshParts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_pointEvaluator(pointEvaluator)
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
{
  m_isInitialized = recvBulk != nullptr;
}

ElementRecvMesh::ElementRecvMesh(stk::mesh::BulkData* recvBulk,
                                 const stk::mesh::FieldBase* coordinateField,
                                 const stk::mesh::EntityRank searchEntityRank,
                                 const stk::mesh::PartVector& recvParts,
                                 const stk::mesh::Selector& activeSelector,
                                 const stk::ParallelMachine recvComm,
                                 std::shared_ptr<PointEvaluatorInterface> pointEvaluator,
                                 double parametricTolerance, double geometricTolerance)
  : m_bulk(recvBulk)
  , m_meta(recvBulk != nullptr ? &recvBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_searchEntityRank(searchEntityRank)
  , m_meshParts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(activeSelector)
  , m_pointEvaluator(pointEvaluator)
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
{
  m_isInitialized = recvBulk != nullptr;
}

void ElementRecvMesh::consistency_check()
{
  STK_ThrowRequireMsg((m_meta != nullptr) && (m_bulk != nullptr),
                       "Meta and/or Bulk Data objects have not been set");

  STK_ThrowRequireMsg(
      (m_searchEntityRank == stk::topology::ELEM_RANK) || (m_searchEntityRank == stk::topology::FACE_RANK) ||
          (m_searchEntityRank == stk::topology::EDGE_RANK && m_meta->spatial_dimension() == 2),
      "Input object type must be {ELEM_RANK} or {FACE_RANK} or {EDGE_RANK with dimension 2}: " << m_searchEntityRank);

  for(const stk::mesh::Part* part : m_meshParts) {
    STK_ThrowRequireMsg(part_has_proper_entity_rank(part),
                    "All destination parts must be {ELEM_RANK} or {FACE_RANK} or {EDGE_RANK with dimension 2}");
  }
}

void ElementRecvMesh::initialize()
{
  m_isInitialized = m_bulk != nullptr;

  consistency_check();
}

void ElementRecvMesh::initialize(stk::mesh::BulkData* recvBulk)
{
  m_bulk = recvBulk;
  m_meta = recvBulk != nullptr ? &recvBulk->mesh_meta_data() : nullptr;
  m_isInitialized = recvBulk != nullptr;

  consistency_check();
}

void ElementRecvMesh::bounding_boxes(std::vector<BoundingBox>& v_box) const
{
  STK_ThrowAssert(m_isInitialized);
  STK_ThrowAssertMsg((nullptr != m_bulk) && (nullptr != m_meta), "Meta and/or Bulk Data objects have not been set");
  const unsigned nDim = m_meta->spatial_dimension();

  stk::mesh::EntityRank rank = get_objects_rank(m_searchEntityRank, m_meshParts);
  stk::mesh::Selector selector = get_objects_selector(*m_meta, m_meshParts, &m_activeSelector);
  stk::mesh::BucketVector const& buckets = m_bulk->get_buckets(rank, selector);

  std::vector<double> coords;

  for(auto&& ib : buckets) {
    stk::mesh::Bucket& b = *ib;

    for(auto entity : b) {

      stk::search::spmd::EntityKeyPair eKey(make_entity_key_pair(m_bulk, entity));

      stk::topology topo = b.topology();
      size_t numPoints = m_pointEvaluator->num_points(eKey, topo);

      for(size_t p = 0; p < numPoints; ++p) {
        m_pointEvaluator->coordinates(eKey, p, coords);

        Point center;
        if(nDim == 2) {
          center = Point(coords[0], coords[1]);
        }
        else {
          center = Point(coords[0], coords[1], coords[2]);
        }

        EntityKey key(eKey, p);
        EntityProc theIdent(key, m_bulk->parallel_rank());
        BoundingBox theBox(Sphere(center, m_searchTolerance), theIdent);

        v_box.push_back(theBox);
      }
    }
  }
  std::sort(v_box.begin(), v_box.end(),
            [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
}

void ElementRecvMesh::coordinates(const EntityKey& k, std::vector<double>& coords) const
{
  STK_ThrowAssert(m_isInitialized);

  size_t p = k.second;

  m_pointEvaluator->coordinates(k.first, p, coords);
}

stk::mesh::EntityId ElementRecvMesh::id(const EntityKey& k) const
{
  STK_ThrowAssert(m_isInitialized);
  return k.first.id();
}

stk::mesh::Entity ElementRecvMesh::entity(const EntityKey& k) const
{
  STK_ThrowAssert(m_isInitialized);
  stk::mesh::Entity e = k.first;
  return e;
}

double ElementRecvMesh::get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& toCoords) const
{
  STK_ThrowAssert(m_isInitialized);
  const stk::mesh::Entity e = k.first;
  return distance_from_nearest_entity_node(*m_bulk, e, m_coordinateField, toCoords);
}

void ElementRecvMesh::centroid(const EntityKey& k, std::vector<double>& centroid) const
{
  STK_ThrowAssert(m_isInitialized);
  stk::mesh::Entity elem = k.first;
  const unsigned nDim = m_coordinateField->mesh_meta_data().spatial_dimension();
  centroid.assign(nDim, 0.0);
  determine_centroid(nDim, elem, *m_coordinateField, centroid.data());
}

void ElementRecvMesh::fill_entity_keys(const stk::mesh::EntityKeyVector& rangeEntities,
                                       std::vector<EntityKey>& elementEntityKeys)
{
  for(stk::mesh::EntityKey key : rangeEntities) {
    stk::mesh::Entity elem = m_bulk->get_entity(key);
    if(!m_bulk->is_valid(elem) || !m_bulk->bucket(elem).owned()) {
      continue;
    }

    stk::search::spmd::EntityKeyPair spmdEntityKey(make_entity_key_pair(m_bulk, elem));

    stk::topology topo = m_bulk->bucket(elem).topology();
    size_t numPoints = m_pointEvaluator->num_points(spmdEntityKey, topo);

    for(size_t p = 0; p < numPoints; ++p) {
      elementEntityKeys.emplace_back(spmdEntityKey, p);
    }
  }
}

} // namespace spmd
} // namespace search
} // namespace stk


