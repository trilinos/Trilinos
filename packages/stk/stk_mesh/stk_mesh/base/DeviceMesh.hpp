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

#ifndef STK_MESH_DEVICEMESH_HPP
#define STK_MESH_DEVICEMESH_HPP

#include <stk_util/stk_config.h>
#include <stk_util/util/StridedArray.hpp>
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include <Kokkos_Core.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <string>
#include <memory>

#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/NgpUtils.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include "stk_mesh/baseImpl/DeviceMeshHostData.hpp"

namespace stk {
namespace mesh {

class DeviceMesh;

struct DeviceBucket {
  using ConnectedNodes    = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedOrdinals = util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
  using Permutations      = util::StridedArray<const stk::mesh::Permutation>;

  KOKKOS_FUNCTION
  DeviceBucket()
    : m_owningMesh(nullptr),
      m_bucketSize(0),
      m_bucketCapacity(0),
      m_bucketId(0),
      m_bucketTopology(),
      m_entityRank(stk::topology::NODE_RANK)
  {}

  KOKKOS_FUNCTION
  unsigned bucket_id() const { return m_bucketId; }

  KOKKOS_FUNCTION
  size_t size() const { return m_bucketSize; }

  KOKKOS_FUNCTION
  size_t capacity() const { return m_bucketCapacity; }

  KOKKOS_FUNCTION
  stk::mesh::EntityRank entity_rank() const { return m_entityRank; }

  KOKKOS_FUNCTION
  stk::topology topology() const { return m_bucketTopology; }

  KOKKOS_INLINE_FUNCTION
  ConnectedEntities get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  KOKKOS_INLINE_FUNCTION
  ConnectedOrdinals get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_edges(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_faces(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_elements(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  stk::mesh::Entity operator[](unsigned offsetIntoBucket) const {
    return m_entities(offsetIntoBucket);
  }

  KOKKOS_FUNCTION
  bool member(stk::mesh::PartOrdinal partOrdinal) const
  {
    for(unsigned i=0; i<m_partOrdinals.size(); i++) {
      if(m_partOrdinals(i) == partOrdinal) {
        return true;
      }
    }
    return false;
  }

  void initialize_bucket_attributes(const stk::mesh::Bucket &bucket);
  void initialize_fixed_data_from_host(const stk::mesh::Bucket &bucket);
  void update_entity_data_from_host(const stk::mesh::Bucket &bucket);

  void resize_device_views(const stk::mesh::Bucket &bucket);
  std::pair<unsigned, unsigned> scan_entities_for_nodal_connectivity(const stk::mesh::Bucket & bucket);

  EntityViewType m_entities;
  BucketConnectivityType m_nodeConnectivity;
  OrdinalViewType m_nodeConnectivityOffsets;

  OrdinalViewType m_nodeOrdinals;

  PartOrdinalViewType m_partOrdinals;

  const stk::mesh::DeviceMesh* m_owningMesh;

  unsigned m_bucketSize;
  unsigned m_bucketCapacity;

  unsigned m_bucketId;
  stk::topology m_bucketTopology;
  stk::mesh::EntityRank m_entityRank;
};

class DeviceMesh : public NgpMeshBase
{
public:
  using MeshExecSpace     = stk::ngp::ExecSpace;
  using ConnectedNodes    = DeviceBucket::ConnectedNodes;
  using ConnectedEntities = DeviceBucket::ConnectedEntities;
  using ConnectedOrdinals = DeviceBucket::ConnectedOrdinals;
  using Permutations      = DeviceBucket::Permutations;
  using MeshIndex         = FastMeshIndex;
  using BucketType        = DeviceBucket;

  KOKKOS_FUNCTION
  DeviceMesh()
    : NgpMeshBase(),
      bulk(nullptr),
      spatial_dimension(0),
      synchronizedCount(0),
      deviceMeshHostData(nullptr)
  {}

  explicit DeviceMesh(const stk::mesh::BulkData& b)
    : NgpMeshBase(),
      bulk(&b),
      spatial_dimension(b.mesh_meta_data().spatial_dimension()),
      synchronizedCount(0),
      endRank(static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count())),
      copyCounter("copy_counter"),
      deviceMeshHostData(nullptr)
  {
    bulk->register_device_mesh();
    deviceMeshHostData = impl::get_ngp_mesh_host_data(*bulk);
    update_mesh();
  }

  KOKKOS_DEFAULTED_FUNCTION DeviceMesh(const DeviceMesh &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMesh(DeviceMesh &&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMesh& operator=(const DeviceMesh &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMesh& operator=(DeviceMesh &&) = default;

  KOKKOS_FUNCTION
  virtual ~DeviceMesh() override {
    clear_buckets_and_views();
  }

  void update_mesh() override;

  KOKKOS_FUNCTION
  unsigned get_spatial_dimension() const
  {
    return spatial_dimension;
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityId identifier(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()].id();
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityRank entity_rank(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()].rank();
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityKey entity_key(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()];
  }

  KOKKOS_FUNCTION
  stk::mesh::Entity get_entity(stk::mesh::EntityRank rank,
                               const stk::mesh::FastMeshIndex& meshIndex) const
  {
    return buckets[rank](meshIndex.bucket_id)[meshIndex.bucket_ord];
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
  {
    if (connectedRank == stk::topology::NODE_RANK)
    {
      return buckets[rank](entity.bucket_id).get_connected_entities(entity.bucket_ord, connectedRank);
    }
  
    int entityOffset = bucketEntityOffsets[rank](entity.bucket_id) + entity.bucket_ord;
    int connectivityOffset = entityConnectivityOffset[rank][connectedRank](entityOffset);
    size_t numConnected = entityConnectivityOffset[rank][connectedRank](entityOffset+1)
        - connectivityOffset;
    ConnectedEntities connectedEntities(nullptr, 0);
    if (numConnected > 0) {
      int stride = 1;
      connectedEntities =
          ConnectedEntities(&(sparseConnectivity[rank][connectedRank](connectivityOffset)), numConnected, stride);
    }
    return connectedEntities;
  }

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return buckets[rank](entity.bucket_id).get_nodes(entity.bucket_ord);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_edges(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_faces(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_elements(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
  {
    if (connectedRank == stk::topology::NODE_RANK) {
      return buckets[rank](entity.bucket_id).get_connected_ordinals(entity.bucket_ord, connectedRank);
    }
  
    int entityOffset = bucketEntityOffsets[rank](entity.bucket_id) + entity.bucket_ord;
    int connectivityOffset = entityConnectivityOffset[rank][connectedRank](entityOffset);
    size_t numConnected = entityConnectivityOffset[rank][connectedRank](entityOffset+1)
        - connectivityOffset;
    ConnectedOrdinals connectedOrdinals(nullptr, 0);
    if (numConnected > 0)
    {
      int stride = 1;
      connectedOrdinals = ConnectedOrdinals(
          &(sparseConnectivityOrdinals[rank][connectedRank](connectivityOffset)), numConnected, stride);
    }
    return connectedOrdinals;
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_node_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_edge_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_face_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_element_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
  {
    Permutations permutations(nullptr, 0);
    if (connectedRank == stk::topology::NODE_RANK)
    {
      return permutations;
    }
  
    int entityOffset = bucketEntityOffsets[rank](entity.bucket_id) + entity.bucket_ord;
    int connectivityOffset = entityConnectivityOffset[rank][connectedRank](entityOffset);
    size_t numConnected = entityConnectivityOffset[rank][connectedRank](entityOffset+1)
        - connectivityOffset;
    if (numConnected > 0)
    {
      int stride = 1;
      permutations = Permutations(&(sparsePermutations[rank][connectedRank](connectivityOffset)), numConnected, stride);
    }
    return permutations;
  }

  KOKKOS_FUNCTION
  Permutations get_node_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_edge_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_face_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_element_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
  {
    return device_mesh_index(entity);
  }

  KOKKOS_FUNCTION
  stk::mesh::FastMeshIndex device_mesh_index(stk::mesh::Entity entity) const
  {
    return deviceMeshIndices(entity.local_offset());
  }

  stk::NgpVector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
  {
    return stk::mesh::get_bucket_ids(get_bulk_on_host(), rank, selector);
  }

  KOKKOS_FUNCTION
  unsigned num_buckets(stk::mesh::EntityRank rank) const
  {
    return buckets[rank].size();
  }

  KOKKOS_FUNCTION
  const DeviceBucket &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
  {
    buckets[rank](index).m_owningMesh = this;
    return buckets[rank](index);
  }

  KOKKOS_FUNCTION
  DeviceCommMapIndices volatile_fast_shared_comm_map(stk::topology::rank_t rank, int proc) const
  {
    const size_t dataBegin = volatileFastSharedCommMapOffset[rank][proc];
    const size_t dataEnd   = volatileFastSharedCommMapOffset[rank][proc+1];
    DeviceCommMapIndices buffer = Kokkos::subview(volatileFastSharedCommMap[rank], Kokkos::pair<size_t, size_t>(dataBegin, dataEnd));
    return buffer;
  }

  void clear()
  {
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<stk::topology::NUM_RANKS; rank++)
      buckets[rank] = BucketView();
  }

  const stk::mesh::BulkData &get_bulk_on_host() const
  {
    STK_ThrowRequireMsg(bulk != nullptr, "DeviceMesh::get_bulk_on_host, bulk==nullptr");
    return *bulk;
  }

  bool is_up_to_date() const {
    if(bulk == nullptr) { return false; }
    return synchronizedCount == bulk->synchronized_count();
  }

private:
  void set_entity_keys(const stk::mesh::BulkData& bulk_in);

  void set_bucket_entity_offsets(const stk::mesh::BulkData& bulk_in);

  void fill_sparse_connectivities(const stk::mesh::BulkData& bulk_in);

  KOKKOS_FUNCTION
  bool is_last_mesh_copy() const
  {
    return (copyCounter.use_count() == 1);
  }

  KOKKOS_FUNCTION
  bool is_last_bucket_reference(unsigned rank = stk::topology::NODE_RANK) const
  {
    return (buckets[rank].use_count() == 1);
  }

  KOKKOS_FUNCTION
  void clear_buckets_and_views()
  {
    KOKKOS_IF_ON_HOST((
      if (is_last_mesh_copy()) {
        bulk->unregister_device_mesh();
      }
    
      if (is_last_bucket_reference()) {
        for (stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++) {
          for (unsigned iBucket = 0; iBucket < buckets[rank].size(); ++iBucket) {
            buckets[rank][iBucket].~DeviceBucket();
          }
        }
      }
    ))
  }

  bool fill_buckets(const stk::mesh::BulkData& bulk_in);

  void fill_mesh_indices(const stk::mesh::BulkData& bulk_in);

  void copy_entity_keys_to_device();

  void copy_mesh_indices_to_device();

  void copy_bucket_entity_offsets_to_device();

  void copy_sparse_connectivities_to_device();

  void copy_volatile_fast_shared_comm_map_to_device();


  using BucketView = Kokkos::View<DeviceBucket*, stk::ngp::UVMMemSpace>;
  const stk::mesh::BulkData *bulk;
  unsigned spatial_dimension;
  unsigned synchronizedCount;
  stk::mesh::EntityRank endRank;
  Kokkos::View<int[1], Kokkos::HostSpace> copyCounter;
  impl::DeviceMeshHostData* deviceMeshHostData;

  EntityKeyViewType entityKeys;

  BucketView buckets[stk::topology::NUM_RANKS];
  HostMeshIndexType hostMeshIndices;
  MeshIndexType deviceMeshIndices;

  BucketEntityOffsetsViewType bucketEntityOffsets[stk::topology::NUM_RANKS];
  UnsignedViewType entityConnectivityOffset[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
  EntityViewType sparseConnectivity[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
  OrdinalViewType sparseConnectivityOrdinals[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
  PermutationViewType sparsePermutations[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
  UnsignedViewType volatileFastSharedCommMapOffset[stk::topology::NUM_RANKS];
  FastSharedCommMapViewType volatileFastSharedCommMap[stk::topology::NUM_RANKS];
};

KOKKOS_INLINE_FUNCTION
DeviceBucket::ConnectedEntities
DeviceBucket::get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  if (connectedRank == stk::topology::NODE_RANK) {
    const unsigned numNodes = m_nodeConnectivityOffsets(offsetIntoBucket+1)-m_nodeConnectivityOffsets(offsetIntoBucket);
    const size_t nodeOffset = m_nodeConnectivityOffsets(offsetIntoBucket);
    return ConnectedEntities(&m_nodeConnectivity(nodeOffset), numNodes, 1);
  }
  STK_NGP_ThrowAssert(m_owningMesh != nullptr);
  stk::mesh::FastMeshIndex meshIndex{bucket_id(), offsetIntoBucket};
  return m_owningMesh->get_connected_entities(entity_rank(), meshIndex, connectedRank);
}

KOKKOS_INLINE_FUNCTION
DeviceBucket::ConnectedOrdinals
DeviceBucket::get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  if (connectedRank == stk::topology::NODE_RANK) {
    const unsigned numNodes = m_nodeConnectivityOffsets(offsetIntoBucket+1)-m_nodeConnectivityOffsets(offsetIntoBucket);
    return ConnectedOrdinals(m_nodeOrdinals.data(), numNodes, 1);
  }
  STK_NGP_ThrowAssert(m_owningMesh != nullptr);
  stk::mesh::FastMeshIndex meshIndex{bucket_id(), offsetIntoBucket};
  return m_owningMesh->get_connected_ordinals(entity_rank(), meshIndex, connectedRank);
}

}
}

#endif

