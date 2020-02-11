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

  void initialize(const stk::mesh::Bucket &bucket) {
    unsigned numNodesPerEntity = bucket.topology().num_nodes();
    bucketId = bucket.bucket_id();
    bucketCapacity = bucket.capacity();
    entityRank = bucket.entity_rank();
    bucketTopology = bucket.topology();
    const std::string bktIdStr(std::to_string(bucketId));

    entities = EntityViewType(Kokkos::ViewAllocateWithoutInitializing("BucketEntities"+bktIdStr), bucket.size());
    hostEntities = Kokkos::create_mirror_view(entities);
    nodeConnectivity = BucketConnectivityType(Kokkos::ViewAllocateWithoutInitializing("BucketConnectivity"+bktIdStr), bucketCapacity, numNodesPerEntity);
    hostNodeConnectivity = Kokkos::create_mirror_view(nodeConnectivity);

    nodeOrdinals = OrdinalViewType(Kokkos::ViewAllocateWithoutInitializing("NodeOrdinals"+bktIdStr), static_cast<size_t>(numNodesPerEntity));
    hostNodeOrdinals = Kokkos::create_mirror_view(nodeOrdinals);
    for(unsigned i=0; i<numNodesPerEntity; ++i)
    {
      hostNodeOrdinals(i) = static_cast<stk::mesh::ConnectivityOrdinal>(i);
    }
    const stk::mesh::PartVector& parts = bucket.supersets();
    partOrdinals = PartOrdinalViewType(Kokkos::ViewAllocateWithoutInitializing("PartOrdinals"+bktIdStr), parts.size());
    hostPartOrdinals = Kokkos::create_mirror_view(partOrdinals);
    for(unsigned i=0; i<parts.size(); ++i)
    {
      hostPartOrdinals(i) = parts[i]->mesh_meta_data_ordinal();
    }
  }

  void copy_to_device()
  {
    Kokkos::deep_copy(entities, hostEntities);
    Kokkos::deep_copy(nodeConnectivity, hostNodeConnectivity);
    Kokkos::deep_copy(nodeOrdinals, hostNodeOrdinals);
    Kokkos::deep_copy(partOrdinals, hostPartOrdinals);
  }

  STK_FUNCTION
  unsigned bucket_id() const { return bucketId; }

  STK_FUNCTION
  size_t size() const { return entities.size(); }

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
  bool member(stk::mesh::PartOrdinal partOrdinal) const {
    for(unsigned i=0; i<partOrdinals.size(); i++) {
      if(partOrdinals(i) == partOrdinal) {
        return true;
      }
    }
    return false;
  }

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
  DeviceMesh() : bulk(nullptr), spatial_dimension(0), synchronizedCount(0) {}

  explicit DeviceMesh(const stk::mesh::BulkData& b)
    : bulk(&b),
      spatial_dimension(b.mesh_meta_data().spatial_dimension()),
      synchronizedCount(b.synchronized_count()),
      endRank(static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count())),
      isInitialized(false)
  {
    update_mesh();
  }

  void update_mesh()
  {
    set_entity_keys(*bulk);
    copy_entity_keys_to_device();
    set_bucket_entity_offsets(*bulk);
    copy_bucket_entity_offsets_to_device();
    fill_sparse_connectivities(*bulk);
    copy_sparse_connectivities_to_device();
    fill_volatile_fast_shared_comm_map(*bulk);
    copy_volatile_fast_shared_comm_map_to_device();

    hostMeshIndices = Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror(Kokkos::ViewAllocateWithoutInitializing("host_mesh_indices"),
                                                                          bulk->get_size_of_entity_index_space());

    if (isInitialized) {
      clear_buckets();
    }
    fill_buckets(*bulk);
    fill_mesh_indices(*bulk);

    copy_mesh_indices_to_device();
    synchronizedCount = bulk->synchronized_count();
    isInitialized = true;
  }

  void update_buckets() const
  {
    //Should this be a throw or a no-op? DeviceMesh certainly can't update buckets,
    //but an app might call this in code that is written for stk::mesh::NgpMesh (i.e., can
    //work for either HostMesh or DeviceMesh.
    //
    //ThrowRequireMsg(false,"ERROR, update_buckets not supported for stk::mesh::DeviceMesh");
  }

  STK_FUNCTION
  DeviceMesh(const DeviceMesh &) = default;

  STK_FUNCTION
  ~DeviceMesh() {
    clear_buckets();
  }

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
  ConnectedNodes get_nodes(const DeviceMeshIndex &entity) const
  {
    return buckets[entity.bucket->entity_rank()](entity.bucket->bucket_id()).get_nodes(entity.bucketOrd);
  }

  STK_FUNCTION
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
    if (numConnected > 0)
    {
      int stride = 1;
      connectedEntities = ConnectedEntities(&sparseConnectivity[rank][connectedRank](connectivityOffset), numConnected, stride);
    }
    return connectedEntities;
  }

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
  ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
  {
    if (connectedRank == stk::topology::NODE_RANK)
    {
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
      connectedOrdinals = ConnectedOrdinals(&sparseConnectivityOrdinals[rank][connectedRank](connectivityOffset), numConnected, stride);
    }
    return connectedOrdinals;
  }

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
      permutations = Permutations(&sparsePermutations[rank][connectedRank](connectivityOffset), numConnected, stride);
    }
    return permutations;
  }

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

  void set_entity_keys(const stk::mesh::BulkData& bulk_in)
  {
    unsigned totalNumEntityKeys = bulk_in.get_size_of_entity_index_space();

    entityKeys = Kokkos::View<stk::mesh::EntityKey*,MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>(Kokkos::ViewAllocateWithoutInitializing("entityKeys"), totalNumEntityKeys);
    hostEntityKeys = Kokkos::create_mirror_view(entityKeys);

    for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
      const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
      for (unsigned i = 0; i < stkBuckets.size(); ++i) {
        const stk::mesh::Bucket & bucket = *stkBuckets[i];
        for (unsigned j = 0; j < bucket.size(); ++j) {
          stk::mesh::Entity entity = bucket[j];
          hostEntityKeys[entity.local_offset()] = bulk_in.entity_key(entity);
        }
      }
    }
  }

  void set_bucket_entity_offsets(const stk::mesh::BulkData& bulk_in)
  {
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
    {
      const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
      bucketEntityOffsets[rank] = Kokkos::View<int*,MemSpace>(Kokkos::ViewAllocateWithoutInitializing("bucketEntityOffsets"), stkBuckets.size()+1);
      hostBucketEntityOffsets[rank] = Kokkos::create_mirror_view(bucketEntityOffsets[rank]);

      int bucketOffsetIntoEntities = 0;
      for(unsigned i=0; i<stkBuckets.size(); ++i)
      {
        hostBucketEntityOffsets[rank](i) = bucketOffsetIntoEntities;
        bucketOffsetIntoEntities += stkBuckets[i]->size();
      }
      hostBucketEntityOffsets[rank](stkBuckets.size()) = bucketOffsetIntoEntities;
    }
  }

  void fill_sparse_connectivities(const stk::mesh::BulkData& bulk_in)
  {
    unsigned totalNumConnectedEntities[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS] = {{0}, {0}, {0}, {0}};
    unsigned totalNumPermutations[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS] = {{0}, {0}, {0}, {0}};
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
    {
      const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
      for(unsigned iBucket=0; iBucket<stkBuckets.size(); ++iBucket)
      {
        const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
        const stk::mesh::EntityRank endConnectedRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
        {
          const bool hasPermutation = stkBucket.has_permutation(connectedRank);
          for(unsigned iEntity=0; iEntity<stkBucket.size(); ++iEntity)
          {
            const unsigned numConnected = stkBucket.num_connectivity(iEntity, connectedRank);
            totalNumConnectedEntities[rank][connectedRank] += numConnected;
            totalNumPermutations[rank][connectedRank] += hasPermutation ? numConnected : 0;
          }
        }
      }

      const stk::mesh::EntityRank endConnectedRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
      for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
      {
        size_t numEntities = hostBucketEntityOffsets[rank](hostBucketEntityOffsets[rank].size()-1);
        entityConnectivityOffset[rank][connectedRank] = UnsignedViewType(Kokkos::ViewAllocateWithoutInitializing("FacesRowMap"), numEntities+1);
        hostEntityConnectivityOffset[rank][connectedRank] = Kokkos::create_mirror_view(entityConnectivityOffset[rank][connectedRank]);

        sparseConnectivity[rank][connectedRank] = EntityViewType(Kokkos::ViewAllocateWithoutInitializing("SparseConnectivity"), totalNumConnectedEntities[rank][connectedRank]);
        hostSparseConnectivity[rank][connectedRank] = Kokkos::create_mirror_view(sparseConnectivity[rank][connectedRank]);

        sparseConnectivityOrdinals[rank][connectedRank] = OrdinalViewType(Kokkos::ViewAllocateWithoutInitializing("SparseConnectivityOrdinals"), totalNumConnectedEntities[rank][connectedRank]);
        hostSparseConnectivityOrdinals[rank][connectedRank] = Kokkos::create_mirror_view(sparseConnectivityOrdinals[rank][connectedRank]);

        sparsePermutations[rank][connectedRank] = PermutationViewType(Kokkos::ViewAllocateWithoutInitializing("SparsePermutations"), totalNumPermutations[rank][connectedRank]);
        hostSparsePermutations[rank][connectedRank] = Kokkos::create_mirror_view(sparsePermutations[rank][connectedRank]);
      }
      int entriesOffsets[stk::topology::NUM_RANKS] = {0};
      unsigned myOffset = 0;
      for(unsigned iBucket=0; iBucket<stkBuckets.size(); ++iBucket)
      {
        const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
        int bucketEntityOffset = hostBucketEntityOffsets[rank](stkBucket.bucket_id());
        for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
        {
          const bool hasPermutation = stkBucket.has_permutation(connectedRank);
          const stk::mesh::Entity* connectedEntities = stkBucket.begin(0, connectedRank);
          const stk::mesh::ConnectivityOrdinal* connectedOrdinals = stkBucket.begin_ordinals(0, connectedRank);
          const stk::mesh::Permutation* permutations = hasPermutation ? stkBucket.begin_permutations(0, connectedRank) : nullptr;
          for(unsigned iEntity=0; iEntity<stkBucket.size(); ++iEntity)
          {
            myOffset = bucketEntityOffset + iEntity;
            unsigned numConnected = stkBucket.num_connectivity(iEntity, connectedRank);

            int entriesOffset = entriesOffsets[connectedRank];
            hostEntityConnectivityOffset[rank][connectedRank](myOffset) = entriesOffset;
            for(unsigned i=0; i<numConnected; ++i)
            {
              hostSparseConnectivity[rank][connectedRank](entriesOffset+i) = connectedEntities[i];
              hostSparseConnectivityOrdinals[rank][connectedRank](entriesOffset+i) = connectedOrdinals[i];
              if (hasPermutation) {
                hostSparsePermutations[rank][connectedRank](entriesOffset+i) = permutations[i];
              }
            }
            entriesOffsets[connectedRank] = entriesOffset + numConnected;
            connectedEntities += numConnected;
            connectedOrdinals += numConnected;
            if (hasPermutation) permutations += numConnected;
          }
        }
      }
      for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
      {
        if (hostEntityConnectivityOffset[rank][connectedRank].size() > myOffset+1)
        {
          hostEntityConnectivityOffset[rank][connectedRank](myOffset+1) = entriesOffsets[connectedRank];
        }
      }
    }
  }

  STK_FUNCTION
  bool is_last_copy() const
  {
    return (buckets[0].use_count() == 1);
  }

  STK_FUNCTION
  void clear_buckets()
  {
    Kokkos::fence();

    if (is_last_copy()) {
      for (stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++) {
        for (unsigned iBucket = 0; iBucket < buckets[rank].size(); ++iBucket) {
          buckets[rank][iBucket].~DeviceBucket();
        }
      }
    }
  }

  void fill_buckets(const stk::mesh::BulkData& bulk_in)
  {
    for (stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++) {
      const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
      unsigned numStkBuckets = stkBuckets.size();
      buckets[rank] = BucketView(Kokkos::ViewAllocateWithoutInitializing("Buckets"), numStkBuckets);

      for (unsigned iBucket = 0; iBucket < numStkBuckets; ++iBucket) {
        const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
        unsigned bucketId = stkBucket.bucket_id();
        unsigned nodesPerElem = stkBucket.topology().num_nodes();

        new (&buckets[rank][bucketId]) DeviceBucket();
        BucketType& bucket = buckets[rank][bucketId];

        bucket.initialize(stkBucket);

        for (unsigned iEntity = 0; iEntity < stkBucket.size(); ++iEntity) {
          stk::mesh::Entity entity = stkBucket[iEntity];
          bucket.hostEntities(iEntity) = entity;
          const stk::mesh::Entity * elemNodes = stkBucket.begin_nodes(iEntity);
          for(unsigned iNode = 0; iNode < nodesPerElem; ++iNode)
          {
            bucket.hostNodeConnectivity(iEntity, iNode) = elemNodes[iNode];
          }
        }
      }

      BucketView &bucketsOfRank = buckets[rank];
      for (unsigned b=0; b<bucketsOfRank.size(); ++b) {
        bucketsOfRank(b).copy_to_device();
      }
    }
  }

  void fill_mesh_indices(const stk::mesh::BulkData& bulk_in)
  {
    for (stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++) {
      const stk::mesh::BucketVector& bkts = bulk_in.buckets(rank);

      for(const stk::mesh::Bucket* bktptr : bkts)
      {
        const stk::mesh::Bucket& bkt = *bktptr;
        const unsigned bktId = bkt.bucket_id();
        for(unsigned i = 0; i < bkt.size(); ++i)
        {
          hostMeshIndices[bkt[i].local_offset()] = stk::mesh::FastMeshIndex{bktId, i};
        }
      }
    }
  }

  void fill_volatile_fast_shared_comm_map(const stk::mesh::BulkData & bulk_in)
  {
    for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::ELEM_RANK; ++rank) {
      if(bulk_in.buckets(rank).size() == 0) { continue; }

      std::vector<size_t> sizePerProc(bulk_in.parallel_size(), 0);

      size_t totalSizeForAllProcs = 0;
      if (bulk_in.parallel_size() > 1) {
        for (int proc = 0; proc < bulk_in.parallel_size(); ++proc) {
          const stk::mesh::BucketIndices & stkBktIndices = bulk_in.volatile_fast_shared_comm_map(rank)[proc];
          sizePerProc[proc] = stkBktIndices.ords.size();
          totalSizeForAllProcs += stkBktIndices.ords.size();
        }
      }

      volatileFastSharedCommMapOffset[rank] = UnsignedViewType(Kokkos::ViewAllocateWithoutInitializing("SharedCommMapOffsets"), sizePerProc.size()+1);
      hostVolatileFastSharedCommMapOffset[rank] = Kokkos::create_mirror_view(volatileFastSharedCommMapOffset[rank]);

      volatileFastSharedCommMap[rank] = FastSharedCommMapViewType(Kokkos::ViewAllocateWithoutInitializing("SharedCommMap"), totalSizeForAllProcs);
      hostVolatileFastSharedCommMap[rank] = Kokkos::create_mirror_view(volatileFastSharedCommMap[rank]);

      size_t entryIndex = 0;
      hostVolatileFastSharedCommMapOffset[rank][0] = 0;
      for (int proc = 0; proc < bulk_in.parallel_size(); ++proc) {
        hostVolatileFastSharedCommMapOffset[rank][proc+1] = hostVolatileFastSharedCommMapOffset[rank][proc] + sizePerProc[proc];

        if (bulk_in.parallel_size() > 1) {
          const stk::mesh::BucketIndices & stkBktIndices = bulk_in.volatile_fast_shared_comm_map(rank)[proc];
          size_t stkOrdinalIndex = 0;
          for (size_t i = 0; i < stkBktIndices.bucket_info.size(); ++i) {
            const unsigned bucketId = stkBktIndices.bucket_info[i].bucket_id;
            const unsigned numEntitiesThisBucket = stkBktIndices.bucket_info[i].num_entities_this_bucket;
            for (size_t n = 0; n < numEntitiesThisBucket; ++n) {
              const unsigned ordinal = stkBktIndices.ords[stkOrdinalIndex++];
              const stk::mesh::FastMeshIndex stkFastMeshIndex{bucketId, ordinal};
              hostVolatileFastSharedCommMap[rank][entryIndex++] = stkFastMeshIndex;
            }
          }
        }
      }
      ThrowRequireMsg(entryIndex == totalSizeForAllProcs, "Unexpected size for volatile fast shared comm map");
    }
  }

  void copy_entity_keys_to_device()
  {
    Kokkos::deep_copy(entityKeys, hostEntityKeys);
  }

  void copy_mesh_indices_to_device()
  {
    unsigned length = hostMeshIndices.size();
    Kokkos::View<stk::mesh::FastMeshIndex*, MemSpace> nonconst_device_mesh_indices(Kokkos::ViewAllocateWithoutInitializing("tmp_dev_mesh_indices"), length);
    Kokkos::deep_copy(nonconst_device_mesh_indices, hostMeshIndices);
    deviceMeshIndices = nonconst_device_mesh_indices;
  }

  void copy_bucket_entity_offsets_to_device()
  {
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
    {
      Kokkos::deep_copy(bucketEntityOffsets[rank], hostBucketEntityOffsets[rank]);
    }
  }

  void copy_sparse_connectivities_to_device()
  {
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
    {
      const stk::mesh::EntityRank endConnectedRank = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
      for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endConnectedRank; connectedRank++)
      {
        Kokkos::deep_copy(entityConnectivityOffset[rank][connectedRank], hostEntityConnectivityOffset[rank][connectedRank]);
        Kokkos::deep_copy(sparseConnectivity[rank][connectedRank], hostSparseConnectivity[rank][connectedRank]);
        Kokkos::deep_copy(sparseConnectivityOrdinals[rank][connectedRank], hostSparseConnectivityOrdinals[rank][connectedRank]);
        Kokkos::deep_copy(sparsePermutations[rank][connectedRank], hostSparsePermutations[rank][connectedRank]);
      }
    }
  }

  void copy_volatile_fast_shared_comm_map_to_device()
  {
    for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::ELEM_RANK; ++rank)
    {
      Kokkos::deep_copy(volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank]);
      Kokkos::deep_copy(volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank]);
    }
  }


  using BucketView = Kokkos::View<DeviceBucket*, UVMMemSpace>;
  const stk::mesh::BulkData *bulk;
  unsigned spatial_dimension;
  unsigned synchronizedCount;
  stk::mesh::EntityRank endRank;
  bool isInitialized;

  EntityKeyViewType::HostMirror hostEntityKeys;
  EntityKeyViewType entityKeys;

  BucketView buckets[stk::topology::NUM_RANKS];
  Kokkos::View<stk::mesh::FastMeshIndex*>::HostMirror hostMeshIndices;
  Kokkos::View<const stk::mesh::FastMeshIndex*, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess> > deviceMeshIndices;

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

