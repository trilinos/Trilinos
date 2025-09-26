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
#include "stk_search_util/spmd/NodeRecvMesh.hpp"
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"
#include "stk_search/DistanceComparison.hpp"
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

NodeRecvMesh::NodeRecvMesh(stk::mesh::BulkData* recvBulk,
                           const stk::mesh::FieldBase* coordinateField,
                           const stk::mesh::PartVector& recvParts,
                           const stk::ParallelMachine recvComm,
                           double parametricTolerance, double geometricTolerance)
  : m_bulk(recvBulk)
  , m_meta(recvBulk != nullptr ? &recvBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_meshParts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
{
  m_isInitialized = recvBulk != nullptr;
}

NodeRecvMesh::NodeRecvMesh(stk::mesh::BulkData* recvBulk,
                           const stk::mesh::FieldBase* coordinateField,
                           const stk::mesh::PartVector& recvParts,
                           const stk::mesh::Selector& activeSelector,
                           const stk::ParallelMachine recvComm,
                           double parametricTolerance, double geometricTolerance)
  : m_bulk(recvBulk)
  , m_meta(recvBulk != nullptr ? &recvBulk->mesh_meta_data() : nullptr)
  , m_coordinateField(coordinateField)
  , m_meshParts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(activeSelector)
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
{
  m_isInitialized = recvBulk != nullptr;
}

void NodeRecvMesh::consistency_check()
{
  STK_ThrowRequireMsg((m_meta != nullptr) && (m_bulk != nullptr),
                       "Meta and/or Bulk Data objects have not been set");

  for(const stk::mesh::Part* part : m_meshParts) {
    STK_ThrowRequireMsg((part != nullptr) && ((part->primary_entity_rank() <= stk::topology::ELEM_RANK) ||
                                              (part->primary_entity_rank() == stk::topology::INVALID_RANK)),
                        "All destination parts must be {ELEM_RANK} or {FACE_RANK} or {EDGE_RANK} or {NODE_RANK}");
  }
}

void NodeRecvMesh::initialize()
{
  m_isInitialized = m_bulk != nullptr;

  consistency_check();
}

void NodeRecvMesh::initialize(stk::mesh::BulkData* recvBulk)
{
  m_bulk = recvBulk;
  m_meta = recvBulk != nullptr ? &recvBulk->mesh_meta_data() : nullptr;
  m_isInitialized = recvBulk != nullptr;

  consistency_check();
}

void NodeRecvMesh::bounding_boxes(std::vector<NodeRecvMesh::BoundingBox>& v_box) const
{
  STK_ThrowAssert(m_isInitialized);
  STK_ThrowAssertMsg((nullptr != m_bulk) && (nullptr != m_meta), "Meta and/or Bulk Data objects have not been set");
  const unsigned nDim = m_meta->spatial_dimension();

  stk::mesh::EntityRank rank = stk::topology::NODE_RANK;
  stk::mesh::Selector selector = get_objects_selector(*m_meta, m_meshParts, &m_activeSelector);
  stk::mesh::BucketVector const& nodeBuckets = m_bulk->get_buckets(rank, selector);

  for(auto&& ib : nodeBuckets) {
    stk::mesh::Bucket& b = *ib;

    for(auto node : b) {
      const double* coords = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, node));

      Point center;
      if(nDim == 2) {
        center = Point(coords[0], coords[1]);
      }
      else {
        center = Point(coords[0], coords[1], coords[2]);
      }

      EntityKey key(make_entity_key_pair(m_bulk,node));
      EntityProc theIdent(key, m_bulk->parallel_rank());
      BoundingBox theBox(Sphere(center, m_searchTolerance), theIdent);
      v_box.push_back(theBox);
    }
  }
  std::sort(v_box.begin(), v_box.end(),
            [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
}

void NodeRecvMesh::coordinates(const EntityKey& k, std::vector<double>& coords) const
{
  STK_ThrowAssert(m_isInitialized);

  stk::mesh::Entity node = k;
  const double* data = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, node));
  const unsigned nDim = m_meta->spatial_dimension();

  coords.assign(data, data+nDim);
}

stk::mesh::EntityId NodeRecvMesh::id(const EntityKey& k) const
{
  return k.id();
}

stk::mesh::Entity NodeRecvMesh::entity(const EntityKey& k) const
{
  STK_ThrowAssert(m_isInitialized);
  stk::mesh::Entity e = k;
  return e;
}

double NodeRecvMesh::get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& toCoords) const
{
  STK_ThrowAssert(m_isInitialized);
  const unsigned nDim = m_meta->spatial_dimension();

  stk::mesh::Entity node = k;
  const double* nodeCoordinates = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, node));

  return stk::search::distance(nDim, nodeCoordinates, toCoords.data());
}

void NodeRecvMesh::centroid(const EntityKey& k, std::vector<double>& centroid) const
{
  stk::mesh::Entity node = k;
  const double* centroidCoord = static_cast<const double *>(stk::mesh::field_data(*m_coordinateField, node));
  const unsigned nDim = m_meta->spatial_dimension();
  centroid.assign(centroidCoord, centroidCoord + nDim);
}

void NodeRecvMesh::fill_entity_keys(const stk::mesh::EntityKeyVector& rangeEntities,
                                    std::vector<EntityKey>& elementEntityKeys)
{
  for(stk::mesh::EntityKey key : rangeEntities) {
    stk::mesh::Entity node = m_bulk->get_entity(key);
    if(!m_bulk->is_valid(node) || !m_bulk->bucket(node).owned()) {
      continue;
    }

    elementEntityKeys.emplace_back(make_entity_key_pair(m_bulk,node));
  }
}

} // namespace spmd
} // namespace search
} // namespace stk
