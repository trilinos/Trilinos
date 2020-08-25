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
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include <Kokkos_Core.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <string>
#include <memory>

#include <stk_mesh/base/NgpSpaces.hpp>
#include <stk_mesh/base/NgpUtils.hpp>
#include <stk_util/util/StkNgpVector.hpp>

namespace stk {
namespace mesh {

using DeviceCommMapIndices = Kokkos::View<stk::mesh::FastMeshIndex*, stk::mesh::MemSpace>;

using EntityKeyViewType         = Kokkos::View<stk::mesh::EntityKey*, MemSpace>;
using EntityViewType            = Kokkos::View<stk::mesh::Entity*, MemSpace>;
using BucketConnectivityType    = Kokkos::View<stk::mesh::Entity**, MemSpace>;
using UnsignedViewType          = Kokkos::View<unsigned*, MemSpace>;
using BoolViewType              = Kokkos::View<bool*, MemSpace>;
using OrdinalViewType           = Kokkos::View<stk::mesh::ConnectivityOrdinal*, MemSpace>;
using PartOrdinalViewType       = Kokkos::View<stk::mesh::PartOrdinal*, MemSpace>;
using PermutationViewType       = Kokkos::View<stk::mesh::Permutation*, MemSpace>;
using FastSharedCommMapViewType = Kokkos::View<stk::mesh::FastMeshIndex*, MemSpace>;
using HostMeshIndexType         = Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror;
using MeshIndexType             = Kokkos::View<const stk::mesh::FastMeshIndex*, MemSpace,
                                               Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

class DeviceMesh;

struct DeviceBucket {
  using ConnectedNodes    = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedOrdinals = util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
  using Permutations      = util::StridedArray<const stk::mesh::Permutation>;

  STK_FUNCTION
  DeviceBucket()
    : bucketId(0), entityRank(stk::topology::NODE_RANK), entities(),
      nodeConnectivity(), hostNodeConnectivity(),
      owningMesh(nullptr), bucketCapacity(0)
  {}

  void initialize_bucket_attributes(const stk::mesh::Bucket &bucket);
  void allocate(const stk::mesh::Bucket &bucket);
  void initialize_from_host(const stk::mesh::Bucket &bucket);
  void update_from_host(const stk::mesh::Bucket &bucket);

  STK_FUNCTION
  unsigned bucket_id() const { return bucketId; }

  STK_FUNCTION
  size_t size() const { return bucketSize; }

  STK_FUNCTION
  stk::mesh::EntityRank entity_rank() const { return entityRank; }

  STK_FUNCTION
  stk::topology topology() const { return bucketTopology; }

  STK_FUNCTION
  unsigned get_num_nodes_per_entity() const { return nodeConnectivity.extent(1); }

  STK_INLINE_FUNCTION
  ConnectedEntities get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  STK_INLINE_FUNCTION
  ConnectedOrdinals get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  STK_FUNCTION
  ConnectedNodes get_nodes(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::NODE_RANK);
  }

  STK_FUNCTION
  ConnectedEntities get_edges(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::EDGE_RANK);
  }

  STK_FUNCTION
  ConnectedEntities get_faces(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::FACE_RANK);
  }

  STK_FUNCTION
  ConnectedEntities get_elements(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::ELEM_RANK);
  }

  STK_FUNCTION
  stk::mesh::Entity operator[](unsigned offsetIntoBucket) const {
    return entities(offsetIntoBucket);
  }

  stk::mesh::Entity host_get_entity(unsigned offsetIntoBucket) const {
    return hostEntities(offsetIntoBucket);
  }

  STK_FUNCTION
  bool member(stk::mesh::PartOrdinal partOrdinal) const;

  unsigned bucketId;
  stk::mesh::EntityRank entityRank;
  stk::topology bucketTopology;

  EntityViewType entities;
  EntityViewType::HostMirror hostEntities;

  BucketConnectivityType nodeConnectivity;
  BucketConnectivityType::HostMirror hostNodeConnectivity;

  OrdinalViewType nodeOrdinals;
  OrdinalViewType::HostMirror hostNodeOrdinals;

  PartOrdinalViewType partOrdinals;
  PartOrdinalViewType::HostMirror hostPartOrdinals;

  const stk::mesh::DeviceMesh* owningMesh;
  unsigned bucketCapacity;
  unsigned bucketSize;
};

struct DeviceMeshIndex
{
  const DeviceBucket *bucket;
  size_t bucketOrd;
};

class DeviceMesh
{
public:
  using MeshExecSpace     = stk::mesh::ExecSpace;
  using ConnectedNodes    = DeviceBucket::ConnectedNodes;
  using ConnectedEntities = DeviceBucket::ConnectedEntities;
  using ConnectedOrdinals = DeviceBucket::ConnectedOrdinals;
  using Permutations      = DeviceBucket::Permutations;
  using MeshIndex         = DeviceMeshIndex;
  using BucketType        = DeviceBucket;

  STK_FUNCTION
  DeviceMesh()
    : bulk(nullptr),
      spatial_dimension(0),
      synchronizedCount(0)
  {}

  explicit DeviceMesh(const stk::mesh::BulkData& b)
    : bulk(&b),
      spatial_dimension(b.mesh_meta_data().spatial_dimension()),
      synchronizedCount(0),
      endRank(static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count())),
      copyCounter("copy_counter")
  {
    bulk->register_device_mesh();
    update_mesh();
  }

  STK_FUNCTION
  DeviceMesh(const DeviceMesh &) = default;

  STK_FUNCTION
  ~DeviceMesh() {
    clear_buckets();
  }

  void update_mesh();

  STK_FUNCTION
  unsigned get_spatial_dimension() const
  {
    return spatial_dimension;
  }

  STK_FUNCTION
  stk::mesh::EntityId identifier(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()].id();
  }

  STK_FUNCTION
  stk::mesh::EntityRank entity_rank(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()].rank();
  }

  STK_FUNCTION
  stk::mesh::EntityKey entity_key(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()];
  }

  STK_FUNCTION
  stk::mesh::Entity get_entity(stk::mesh::EntityRank rank,
                               const stk::mesh::FastMeshIndex& meshIndex) const
  {
    return buckets[rank](meshIndex.bucket_id)[meshIndex.bucket_ord];
  }

  STK_FUNCTION
  ConnectedNodes get_nodes(const DeviceMeshIndex &entity) const
  {
    return buckets[entity.bucket->entity_rank()](entity.bucket->bucket_id()).get_nodes(entity.bucketOrd);
  }

  STK_FUNCTION
  ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const;

  STK_FUNCTION
  ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return buckets[rank](entity.bucket_id).get_nodes(entity.bucket_ord);
  }

  STK_FUNCTION
  ConnectedEntities get_edges(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::EDGE_RANK);
  }

  STK_FUNCTION
  ConnectedEntities get_faces(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::FACE_RANK);
  }

  STK_FUNCTION
  ConnectedEntities get_elements(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::ELEM_RANK);
  }

  STK_FUNCTION
  ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const;

  STK_FUNCTION
  ConnectedOrdinals get_node_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::NODE_RANK);
  }

  STK_FUNCTION
  ConnectedOrdinals get_edge_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::EDGE_RANK);
  }

  STK_FUNCTION
  ConnectedOrdinals get_face_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::FACE_RANK);
  }

  STK_FUNCTION
  ConnectedOrdinals get_element_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::ELEM_RANK);
  }

  STK_FUNCTION
  Permutations get_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const;

  STK_FUNCTION
  Permutations get_node_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::NODE_RANK);
  }

  STK_FUNCTION
  Permutations get_edge_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::EDGE_RANK);
  }

  STK_FUNCTION
  Permutations get_face_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::FACE_RANK);
  }

  STK_FUNCTION
  Permutations get_element_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::ELEM_RANK);
  }

  STK_FUNCTION
  stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
  {
    return device_mesh_index(entity);
  }

  STK_FUNCTION
  stk::mesh::FastMeshIndex device_mesh_index(stk::mesh::Entity entity) const
  {
    return deviceMeshIndices(entity.local_offset());
  }

  const stk::mesh::FastMeshIndex& host_mesh_index(stk::mesh::Entity entity) const
  {
    return hostMeshIndices(entity.local_offset());
  }

  stk::NgpVector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
  {
    return stk::mesh::get_bucket_ids(get_bulk_on_host(), rank, selector);
  }

  STK_FUNCTION
  unsigned num_buckets(stk::mesh::EntityRank rank) const
  {
    return buckets[rank].size();
  }

  STK_FUNCTION
  const DeviceBucket &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
  {
    buckets[rank](index).owningMesh = this;
    return buckets[rank](index);
  }

  STK_FUNCTION
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

  STK_FUNCTION
  bool is_last_mesh_copy() const
  {
    return (copyCounter.use_count() == 1);
  }

  STK_FUNCTION
  bool is_last_bucket_reference(unsigned rank = stk::topology::NODE_RANK) const
  {
    return (buckets[rank].use_count() == 1);
  }

  STK_FUNCTION
  void clear_buckets();

  void fill_buckets(const stk::mesh::BulkData& bulk_in);

  void fill_mesh_indices(const stk::mesh::BulkData& bulk_in);

  void fill_volatile_fast_shared_comm_map(const stk::mesh::BulkData & bulk_in);

  void copy_entity_keys_to_device();

  void copy_mesh_indices_to_device();

  void copy_bucket_entity_offsets_to_device();

  void copy_sparse_connectivities_to_device();

  void copy_volatile_fast_shared_comm_map_to_device();


  using BucketView = Kokkos::View<DeviceBucket*, UVMMemSpace>;
  const stk::mesh::BulkData *bulk;
  unsigned spatial_dimension;
  unsigned synchronizedCount;
  stk::mesh::EntityRank endRank;
  Kokkos::View<int[1], Kokkos::HostSpace> copyCounter;

  EntityKeyViewType::HostMirror hostEntityKeys;
  EntityKeyViewType entityKeys;

  BucketView buckets[stk::topology::NUM_RANKS];
  HostMeshIndexType hostMeshIndices;
  MeshIndexType deviceMeshIndices;

  Kokkos::View<int*,MemSpace> bucketEntityOffsets[stk::topology::NUM_RANKS];
  Kokkos::View<int*,MemSpace>::HostMirror hostBucketEntityOffsets[stk::topology::NUM_RANKS];

  UnsignedViewType entityConnectivityOffset[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
  UnsignedViewType::HostMirror hostEntityConnectivityOffset[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];

  EntityViewType sparseConnectivity[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
  EntityViewType::HostMirror hostSparseConnectivity[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];

  OrdinalViewType sparseConnectivityOrdinals[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
  OrdinalViewType::HostMirror hostSparseConnectivityOrdinals[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];

  PermutationViewType sparsePermutations[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];
  PermutationViewType::HostMirror hostSparsePermutations[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS];

  UnsignedViewType volatileFastSharedCommMapOffset[stk::topology::NUM_RANKS];
  UnsignedViewType::HostMirror hostVolatileFastSharedCommMapOffset[stk::topology::NUM_RANKS];

  FastSharedCommMapViewType volatileFastSharedCommMap[stk::topology::NUM_RANKS];
  FastSharedCommMapViewType::HostMirror hostVolatileFastSharedCommMap[stk::topology::NUM_RANKS];
};

STK_INLINE_FUNCTION
DeviceBucket::ConnectedEntities
DeviceBucket::get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  if (connectedRank == stk::topology::NODE_RANK) {
    return ConnectedEntities(&nodeConnectivity(offsetIntoBucket,0), nodeConnectivity.extent(1), bucketCapacity);
  }
  NGP_ThrowAssert(owningMesh != nullptr);
  stk::mesh::FastMeshIndex meshIndex{bucket_id(), offsetIntoBucket};
  return owningMesh->get_connected_entities(entity_rank(), meshIndex, connectedRank);
}

STK_INLINE_FUNCTION
DeviceBucket::ConnectedOrdinals
DeviceBucket::get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  if (connectedRank == stk::topology::NODE_RANK) {
    return ConnectedOrdinals(&nodeOrdinals(0), nodeOrdinals.extent(0), bucketCapacity);
  }
  NGP_ThrowAssert(owningMesh != nullptr);
  stk::mesh::FastMeshIndex meshIndex{bucket_id(), offsetIntoBucket};
  return owningMesh->get_connected_ordinals(entity_rank(), meshIndex, connectedRank);
}

}
}

#endif

