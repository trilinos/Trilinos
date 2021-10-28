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

#include "stk_mesh/base/DeviceMesh.hpp"

namespace stk {
namespace mesh {

void DeviceBucket::initialize_bucket_attributes(const stk::mesh::Bucket &bucket)
{
  bucketId = bucket.bucket_id();
  bucketCapacity = bucket.capacity();
  bucketSize = bucket.size();
  entityRank = bucket.entity_rank();
  bucketTopology = bucket.topology();
}

void DeviceBucket::allocate(const stk::mesh::Bucket &bucket)
{
  nodeOffsets = OrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "NodeOffsets"), bucket.size()+1);
  hostNodeOffsets = Kokkos::create_mirror_view(Kokkos::HostSpace(), nodeOffsets, Kokkos::WithoutInitializing);

  unsigned maxNodesPerEntity = bucketTopology.num_nodes();
  unsigned totalNumNodes = bucketTopology.num_nodes()*bucketCapacity;
  if (bucket.topology() == stk::topology::INVALID_TOPOLOGY) {
    maxNodesPerEntity = 0;
    totalNumNodes = 0;
    for(unsigned i=0; i<bucket.size(); ++i) {
      totalNumNodes += bucket.num_nodes(i);
      maxNodesPerEntity = std::max(maxNodesPerEntity, bucket.num_nodes(i));
    }
  }

  const stk::mesh::PartVector& parts = bucket.supersets();

  entities = EntityViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "BucketEntities"), bucketCapacity);
  hostEntities = Kokkos::create_mirror_view(Kokkos::HostSpace(), entities, Kokkos::WithoutInitializing);

  nodeConnectivity = BucketConnectivityType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "BucketConnectivity"), totalNumNodes);
  hostNodeConnectivity = Kokkos::create_mirror_view(Kokkos::HostSpace(), nodeConnectivity, Kokkos::WithoutInitializing);

  nodeOrdinals = OrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "NodeOrdinals"),
                                 static_cast<size_t>(maxNodesPerEntity));
  hostNodeOrdinals = Kokkos::create_mirror_view(Kokkos::HostSpace(), nodeOrdinals, Kokkos::WithoutInitializing);

  partOrdinals = PartOrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PartOrdinals"), parts.size());
  hostPartOrdinals = Kokkos::create_mirror_view(Kokkos::HostSpace(), partOrdinals, Kokkos::WithoutInitializing);
}

void DeviceBucket::initialize_from_host(const stk::mesh::Bucket &bucket)
{
  unsigned maxNodesPerEntity = bucket.topology().num_nodes();
  if (bucket.topology() == stk::topology::INVALID_TOPOLOGY) {
    maxNodesPerEntity = 0;
    for(unsigned i=0; i<bucket.size(); ++i) {
      maxNodesPerEntity = std::max(maxNodesPerEntity, bucket.num_nodes(i));
    }
  }

  for (unsigned i = 0; i < maxNodesPerEntity; ++i) {
    hostNodeOrdinals(i) = static_cast<stk::mesh::ConnectivityOrdinal>(i);
  }
  Kokkos::deep_copy(nodeOrdinals, hostNodeOrdinals);

  const stk::mesh::PartVector& parts = bucket.supersets();
  for (unsigned i = 0; i < parts.size(); ++i) {
    hostPartOrdinals(i) = parts[i]->mesh_meta_data_ordinal();
  }
  Kokkos::deep_copy(partOrdinals, hostPartOrdinals);

  update_from_host(bucket);
}

void DeviceBucket::update_from_host(const stk::mesh::Bucket &bucket)
{
  bucketSize = bucket.size();

  if (bucketSize+1 != hostNodeOffsets.size()) {
    nodeOffsets = OrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "NodeOffsets"), bucketSize+1);
    hostNodeOffsets = Kokkos::create_mirror_view(Kokkos::HostSpace(), nodeOffsets, Kokkos::WithoutInitializing);
  }

  unsigned totalNumNodes = bucket.topology().num_nodes()*bucketCapacity;
  unsigned maxNodesPerEntity = bucket.topology().num_nodes();
  if (bucket.topology() == stk::topology::INVALID_TOPOLOGY) {
    totalNumNodes = 0;
    maxNodesPerEntity = 0;
    for(unsigned i=0; i<bucket.size(); ++i) {
      totalNumNodes += bucket.num_nodes(i);
      maxNodesPerEntity = std::max(maxNodesPerEntity, bucket.num_nodes(i));
    }
  }

  if (totalNumNodes != hostNodeConnectivity.size()) {
    nodeConnectivity = BucketConnectivityType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "BucketConnectivity"), totalNumNodes);
    hostNodeConnectivity = Kokkos::create_mirror_view(Kokkos::HostSpace(), nodeConnectivity, Kokkos::WithoutInitializing);
  }

  if (maxNodesPerEntity != hostNodeOrdinals.size()) {
    nodeOrdinals = OrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "NodeOrdinals"),
                                   static_cast<size_t>(maxNodesPerEntity));
    hostNodeOrdinals = Kokkos::create_mirror_view(Kokkos::HostSpace(), nodeOrdinals, Kokkos::WithoutInitializing);
    for (unsigned i = 0; i < maxNodesPerEntity; ++i) {
      hostNodeOrdinals(i) = static_cast<stk::mesh::ConnectivityOrdinal>(i);
    }
    Kokkos::deep_copy(nodeOrdinals, hostNodeOrdinals);
  }

  unsigned nodeOffset = 0;
  for (unsigned iEntity = 0; iEntity < bucket.size(); ++iEntity) {
    stk::mesh::Entity entity = bucket[iEntity];
    hostEntities(iEntity) = entity;
    const unsigned nodesPerEntity = bucket.num_nodes(iEntity);
    const stk::mesh::Entity * elemNodes = bucket.begin_nodes(iEntity);
    for (unsigned iNode = 0; iNode < nodesPerEntity; ++iNode) {
      hostNodeConnectivity(nodeOffset + iNode) = elemNodes[iNode];
    }
    hostNodeOffsets(iEntity) = nodeOffset;
    nodeOffset += nodesPerEntity;
  }
  hostNodeOffsets(bucket.size()) = nodeOffset;

  Kokkos::deep_copy(entities, hostEntities);
  Kokkos::deep_copy(nodeConnectivity, hostNodeConnectivity);
  Kokkos::deep_copy(nodeOffsets, hostNodeOffsets);
}

void DeviceMesh::update_mesh()
{
  if (is_up_to_date()) return;

  require_ngp_mesh_rank_limit(bulk->mesh_meta_data());

  const bool anyChanges = fill_buckets(*bulk);

  if (anyChanges) {
    set_entity_keys(*bulk);
    copy_entity_keys_to_device();
    set_bucket_entity_offsets(*bulk);
    copy_bucket_entity_offsets_to_device();
    fill_sparse_connectivities(*bulk);
    copy_sparse_connectivities_to_device();
    fill_volatile_fast_shared_comm_map(*bulk);
    copy_volatile_fast_shared_comm_map_to_device();
    fill_mesh_indices(*bulk);
    copy_mesh_indices_to_device();
  }

  synchronizedCount = bulk->synchronized_count();
}

bool DeviceMesh::fill_buckets(const stk::mesh::BulkData& bulk_in)
{
  bool anyBucketChanges = false;

  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
    const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
    unsigned numStkBuckets = stkBuckets.size();

    BucketView bucketBuffer(Kokkos::view_alloc(Kokkos::WithoutInitializing, "BucketBuffer"), numStkBuckets);

    if (numStkBuckets != buckets[rank].size()) {
      anyBucketChanges = true;
    }

    for (unsigned iBucket = 0; iBucket < numStkBuckets; ++iBucket) {
      stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
      const unsigned ngpBucketId = stkBucket.ngp_bucket_id();

      if (ngpBucketId == INVALID_BUCKET_ID) {
        // New bucket on host
        new (&bucketBuffer[iBucket]) DeviceBucket();
        bucketBuffer[iBucket].initialize_bucket_attributes(stkBucket);
        bucketBuffer[iBucket].allocate(stkBucket);
        bucketBuffer[iBucket].initialize_from_host(stkBucket);
        anyBucketChanges = true;
      }
      else {
        // Pre-existing bucket on host
        new (&bucketBuffer[iBucket]) DeviceBucket(buckets[rank][ngpBucketId]);
        if (stkBucket.is_modified()) {
          bucketBuffer[iBucket].update_from_host(stkBucket);
          anyBucketChanges = true;
        }
        bucketBuffer[iBucket].bucketId = stkBucket.bucket_id();
      }

      stkBucket.set_ngp_bucket_id(iBucket);
    }

    if (is_last_bucket_reference(rank)) {
      for (unsigned iBucket = 0; iBucket < buckets[rank].size(); ++iBucket) {
        buckets[rank][iBucket].~DeviceBucket();
      }
    }

    buckets[rank] = bucketBuffer;
  }

  return anyBucketChanges;
}

constexpr double RESIZE_FACTOR = 0.05;

template <typename DEVICE_VIEW, typename HOST_VIEW>
inline void reallocate_views(DEVICE_VIEW & deviceView, HOST_VIEW & hostView, size_t requiredSize, double resizeFactor = 0.0)
{
  const size_t currentSize = deviceView.extent(0);
  const size_t shrinkThreshold = currentSize - static_cast<size_t>(2*resizeFactor*currentSize);
  const bool needGrowth = (requiredSize > currentSize);
  const bool needShrink = (requiredSize < shrinkThreshold);

  if (needGrowth || needShrink) {
    const size_t newSize = requiredSize + static_cast<size_t>(resizeFactor*requiredSize);
    deviceView = DEVICE_VIEW(Kokkos::view_alloc(Kokkos::WithoutInitializing, deviceView.label()), newSize);
    hostView = Kokkos::create_mirror_view(Kokkos::HostSpace(), deviceView, Kokkos::WithoutInitializing);
  }
}

void DeviceMesh::set_entity_keys(const stk::mesh::BulkData& bulk_in)
{
  unsigned totalNumEntityKeys = bulk_in.get_size_of_entity_index_space();

  reallocate_views(entityKeys, hostEntityKeys, totalNumEntityKeys, RESIZE_FACTOR);

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

void DeviceMesh::set_bucket_entity_offsets(const stk::mesh::BulkData& bulk_in)
{
  for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
  {
    const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
    reallocate_views(bucketEntityOffsets[rank], hostBucketEntityOffsets[rank], stkBuckets.size()+1, RESIZE_FACTOR);

    int bucketOffsetIntoEntities = 0;
    for (unsigned i = 0; i < stkBuckets.size(); ++i)
    {
      hostBucketEntityOffsets[rank](i) = bucketOffsetIntoEntities;
      bucketOffsetIntoEntities += stkBuckets[i]->size();
    }
    for (unsigned i = stkBuckets.size(); i < hostBucketEntityOffsets[rank].extent(0); ++i) {
      hostBucketEntityOffsets[rank](i) = bucketOffsetIntoEntities;
    }
  }
}

void DeviceMesh::fill_sparse_connectivities(const stk::mesh::BulkData& bulk_in)
{
  unsigned totalNumConnectedEntities[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS] = {{0}, {0}, {0}, {0}, {0}};
  unsigned totalNumPermutations[stk::topology::NUM_RANKS][stk::topology::NUM_RANKS] = {{0}, {0}, {0}, {0}, {0}};
  for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
  {
    const stk::mesh::BucketVector& stkBuckets = bulk_in.buckets(rank);
    for(unsigned iBucket=0; iBucket<stkBuckets.size(); ++iBucket)
    {
      const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
      for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; connectedRank++)
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

    for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; connectedRank++)
    {
      size_t numEntities = hostBucketEntityOffsets[rank](hostBucketEntityOffsets[rank].size()-1);
      reallocate_views(entityConnectivityOffset[rank][connectedRank], hostEntityConnectivityOffset[rank][connectedRank],
                       numEntities+1, RESIZE_FACTOR);

      reallocate_views(sparseConnectivity[rank][connectedRank], hostSparseConnectivity[rank][connectedRank],
                       totalNumConnectedEntities[rank][connectedRank], RESIZE_FACTOR);

      reallocate_views(sparseConnectivityOrdinals[rank][connectedRank], hostSparseConnectivityOrdinals[rank][connectedRank],
                       totalNumConnectedEntities[rank][connectedRank], RESIZE_FACTOR);

      reallocate_views(sparsePermutations[rank][connectedRank], hostSparsePermutations[rank][connectedRank],
                       totalNumPermutations[rank][connectedRank], RESIZE_FACTOR);
    }
    int entriesOffsets[stk::topology::NUM_RANKS] = {0};
    unsigned myOffset = 0;
    for(unsigned iBucket=0; iBucket<stkBuckets.size(); ++iBucket)
    {
      const stk::mesh::Bucket& stkBucket = *stkBuckets[iBucket];
      int bucketEntityOffset = hostBucketEntityOffsets[rank](stkBucket.bucket_id());
      for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; connectedRank++)
      {
        const bool hasPermutation = stkBucket.has_permutation(connectedRank);
        for(unsigned iEntity=0; iEntity<stkBucket.size(); ++iEntity)
        {
          myOffset = bucketEntityOffset + iEntity;
          unsigned numConnected = stkBucket.num_connectivity(iEntity, connectedRank);
          const stk::mesh::Entity* connectedEntities = stkBucket.begin(iEntity, connectedRank);
          const stk::mesh::ConnectivityOrdinal* connectedOrdinals = stkBucket.begin_ordinals(iEntity, connectedRank);
          const stk::mesh::Permutation* permutations = hasPermutation ? stkBucket.begin_permutations(iEntity, connectedRank) : nullptr;

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
        }
      }
    }
    for (stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; connectedRank++)
    {
      for (unsigned i = myOffset+1; i < hostEntityConnectivityOffset[rank][connectedRank].extent(0); ++i) {
        hostEntityConnectivityOffset[rank][connectedRank](i) = entriesOffsets[connectedRank];
      }
    }
  }
}

void DeviceMesh::fill_mesh_indices(const stk::mesh::BulkData& bulk_in)
{
  const size_t indexSpaceSize = bulk->get_size_of_entity_index_space();
  hostMeshIndices = HostMeshIndexType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "host_mesh_indices"), indexSpaceSize);

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

void DeviceMesh::fill_volatile_fast_shared_comm_map(const stk::mesh::BulkData & bulk_in)
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

    reallocate_views(volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank],
                     sizePerProc.size()+1, RESIZE_FACTOR);

    reallocate_views(volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank],
                     totalSizeForAllProcs, RESIZE_FACTOR);

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

void DeviceMesh::copy_entity_keys_to_device()
{
  Kokkos::deep_copy(entityKeys, hostEntityKeys);
}

void DeviceMesh::copy_mesh_indices_to_device()
{
  unsigned length = hostMeshIndices.size();
  Kokkos::View<stk::mesh::FastMeshIndex*, stk::ngp::MemSpace> nonconst_device_mesh_indices(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_dev_mesh_indices"), length);
  Kokkos::deep_copy(nonconst_device_mesh_indices, hostMeshIndices);
  deviceMeshIndices = nonconst_device_mesh_indices;
}

void DeviceMesh::copy_bucket_entity_offsets_to_device()
{
  for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
  {
    Kokkos::deep_copy(bucketEntityOffsets[rank], hostBucketEntityOffsets[rank]);
  }
}

void DeviceMesh::copy_sparse_connectivities_to_device()
{
  for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
  {
    for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; connectedRank++)
    {
      Kokkos::deep_copy(entityConnectivityOffset[rank][connectedRank], hostEntityConnectivityOffset[rank][connectedRank]);
      Kokkos::deep_copy(sparseConnectivity[rank][connectedRank], hostSparseConnectivity[rank][connectedRank]);
      Kokkos::deep_copy(sparseConnectivityOrdinals[rank][connectedRank], hostSparseConnectivityOrdinals[rank][connectedRank]);
      Kokkos::deep_copy(sparsePermutations[rank][connectedRank], hostSparsePermutations[rank][connectedRank]);
    }
  }
}

void DeviceMesh::copy_volatile_fast_shared_comm_map_to_device()
{
  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::ELEM_RANK; ++rank)
  {
    Kokkos::deep_copy(volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank]);
    Kokkos::deep_copy(volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank]);
  }
}

}
}

