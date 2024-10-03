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
  m_bucketId = bucket.bucket_id();
  m_bucketCapacity = bucket.capacity();
  m_bucketSize = bucket.size();
  m_entityRank = bucket.entity_rank();
  m_bucketTopology = bucket.topology();
}

void DeviceBucket::initialize_fixed_data_from_host(const stk::mesh::Bucket &bucket)
{
  const stk::mesh::PartVector& parts = bucket.supersets();
  m_partOrdinals = PartOrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PartOrdinals"),
                                       parts.size());
  auto hostPartOrdinals = HostPartOrdinalViewType(bucket.superset_part_ordinals().first, parts.size());
  Kokkos::deep_copy(m_partOrdinals, hostPartOrdinals);
}

std::pair<unsigned, unsigned>
DeviceBucket::scan_entities_for_nodal_connectivity(const stk::mesh::Bucket & bucket)
{
  if (bucket.topology() == stk::topology::INVALID_TOPOLOGY) {
    unsigned maxNodesPerEntity = 0;
    unsigned totalNumConnectedNodes = 0;
    for (unsigned i = 0; i < bucket.size(); ++i) {
      maxNodesPerEntity = std::max(maxNodesPerEntity, bucket.num_nodes(i));
      totalNumConnectedNodes += bucket.num_nodes(i);
    }
    return std::make_pair(maxNodesPerEntity, totalNumConnectedNodes);
  }

  return std::make_pair<unsigned, unsigned>(bucket.topology().num_nodes(),
                                            bucket.topology().num_nodes() * m_bucketCapacity);
}

void DeviceBucket::resize_device_views(const stk::mesh::Bucket & bucket)
{
  Kokkos::Profiling::pushRegion("resize_device_views()");

  const auto [maxNodesPerEntity, totalNumConnectedNodes] = scan_entities_for_nodal_connectivity(bucket);

  if (m_nodeOrdinals.size() != maxNodesPerEntity) {
    m_nodeOrdinals = OrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "NodeOrdinals"),
                                     static_cast<size_t>(maxNodesPerEntity));
    OrdinalViewType& nodeOrds = m_nodeOrdinals; //local var to avoid implicit this capture
    Kokkos::parallel_for(Kokkos::RangePolicy<stk::ngp::ExecSpace>(0, maxNodesPerEntity),
      KOKKOS_LAMBDA(const int i) {
        nodeOrds(i) = static_cast<stk::mesh::ConnectivityOrdinal>(i);
      });
  }

  if (m_entities.size() != m_bucketCapacity) {
    m_entities = EntityViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "BucketEntities"), m_bucketCapacity);
    STK_ThrowRequireMsg(m_bucketCapacity > 0, "bucket capacity must be greater than 0");
  }

  if (m_nodeConnectivity.size() != totalNumConnectedNodes) {
    m_nodeConnectivity = BucketConnectivityType(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                                   "NodeConnectivity"), totalNumConnectedNodes);
  }

  if (m_nodeConnectivityOffsets.size() != m_bucketCapacity+1) {
    m_nodeConnectivityOffsets = OrdinalViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                                   "NodeConnectivityOffsets"), m_bucketCapacity+1);
  }
  Kokkos::Profiling::popRegion();
}

void DeviceBucket::update_entity_data_from_host(const stk::mesh::Bucket &bucket)
{
  Kokkos::Profiling::pushRegion("update_entity_data_from_host()");

  m_bucketSize = bucket.size();
  m_bucketCapacity = bucket.capacity();

  resize_device_views(bucket);

  Kokkos::Profiling::pushRegion("filling host-side Views");
  auto hostEntities = HostEntityViewType(bucket.begin(), m_bucketCapacity);
  auto hostNodeConnectivity = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeConnectivity);
  auto hostNodeConnectivityOffsets = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeConnectivityOffsets);
  unsigned nodeOffset = 0;
  for (unsigned iEntity = 0; iEntity < bucket.size(); ++iEntity) {
    const unsigned nodesPerEntity = bucket.num_nodes(iEntity);
    const stk::mesh::Entity * elemNodes = bucket.begin_nodes(iEntity);
    for (unsigned iNode = 0; iNode < nodesPerEntity; ++iNode) {
      hostNodeConnectivity(nodeOffset + iNode) = elemNodes[iNode];
    }
    hostNodeConnectivityOffsets(iEntity) = nodeOffset;
    nodeOffset += nodesPerEntity;
  }
  hostNodeConnectivityOffsets(bucket.size()) = nodeOffset;
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("deep_copy entities/connectivity/offsets");
  Kokkos::deep_copy(m_entities, hostEntities);
  Kokkos::deep_copy(m_nodeConnectivity, hostNodeConnectivity);
  Kokkos::deep_copy(m_nodeConnectivityOffsets, hostNodeConnectivityOffsets);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::popRegion();
}

void DeviceMesh::update_mesh()
{
  if (is_up_to_date()) {
    return;
  }

  require_ngp_mesh_rank_limit(bulk->mesh_meta_data());

  Kokkos::Profiling::pushRegion("DeviceMesh::update_mesh");
  const bool anyChanges = fill_buckets(*bulk);

  if (anyChanges) {
    set_entity_keys(*bulk);
    copy_entity_keys_to_device();
    set_bucket_entity_offsets(*bulk);
    copy_bucket_entity_offsets_to_device();
    fill_sparse_connectivities(*bulk);
    copy_sparse_connectivities_to_device();
    copy_volatile_fast_shared_comm_map_to_device();
    fill_mesh_indices(*bulk);
    copy_mesh_indices_to_device();
  }

  synchronizedCount = bulk->synchronized_count();
  Kokkos::Profiling::popRegion();
}

bool DeviceMesh::fill_buckets(const stk::mesh::BulkData& bulk_in)
{
  bool anyBucketChanges = false;

  Kokkos::Profiling::pushRegion("fill_buckets");
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
        Kokkos::Profiling::pushRegion("new bucket");
        // New bucket on host
        new (&bucketBuffer[iBucket]) DeviceBucket();
        bucketBuffer[iBucket].initialize_bucket_attributes(stkBucket);
        bucketBuffer[iBucket].initialize_fixed_data_from_host(stkBucket);
        bucketBuffer[iBucket].update_entity_data_from_host(stkBucket);
        anyBucketChanges = true;
        Kokkos::Profiling::popRegion();
      }
      else {
        Kokkos::Profiling::pushRegion("pre-existing bucket");
        // Pre-existing bucket on host
        new (&bucketBuffer[iBucket]) DeviceBucket(buckets[rank][ngpBucketId]);
        if (stkBucket.is_modified()) {
          bucketBuffer[iBucket].update_entity_data_from_host(stkBucket);
          anyBucketChanges = true;
        }
        bucketBuffer[iBucket].m_bucketId = stkBucket.bucket_id();
        Kokkos::Profiling::popRegion();
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
  Kokkos::Profiling::popRegion();

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
    hostView = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, deviceView);
  }
}

void DeviceMesh::set_entity_keys(const stk::mesh::BulkData& bulk_in)
{
  unsigned totalNumEntityKeys = bulk_in.get_size_of_entity_index_space();
  auto& hostEntityKeys = deviceMeshHostData->hostEntityKeys;

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
  auto& hostBucketEntityOffsets = deviceMeshHostData->hostBucketEntityOffsets;

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
  auto& hostEntityConnectivityOffset = deviceMeshHostData->hostEntityConnectivityOffset;
  auto& hostBucketEntityOffsets = deviceMeshHostData->hostBucketEntityOffsets;
  auto& hostSparseConnectivity = deviceMeshHostData->hostSparseConnectivity;
  auto& hostSparseConnectivityOrdinals = deviceMeshHostData->hostSparseConnectivityOrdinals;
  auto& hostSparsePermutations = deviceMeshHostData->hostSparsePermutations;

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

          int entriesOffset = entriesOffsets[connectedRank];
          hostEntityConnectivityOffset[rank][connectedRank](myOffset) = entriesOffset;

          if (numConnected > 0) {

            const stk::mesh::Entity* connectedEntities = stkBucket.begin(iEntity, connectedRank);
            const stk::mesh::ConnectivityOrdinal* connectedOrdinals = stkBucket.begin_ordinals(iEntity, connectedRank);
            const stk::mesh::Permutation* permutations = hasPermutation ? stkBucket.begin_permutations(iEntity, connectedRank) : nullptr;
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

void DeviceMesh::copy_entity_keys_to_device()
{
  auto& hostEntityKeys = deviceMeshHostData->hostEntityKeys;

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
  auto& hostBucketEntityOffsets = deviceMeshHostData->hostBucketEntityOffsets;

  for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; rank++)
  {
    Kokkos::deep_copy(bucketEntityOffsets[rank], hostBucketEntityOffsets[rank]);
  }
}

void DeviceMesh::copy_sparse_connectivities_to_device()
{
  auto& hostEntityConnectivityOffset = deviceMeshHostData->hostEntityConnectivityOffset;
  auto& hostSparseConnectivity = deviceMeshHostData->hostSparseConnectivity;
  auto& hostSparseConnectivityOrdinals = deviceMeshHostData->hostSparseConnectivityOrdinals;
  auto& hostSparsePermutations = deviceMeshHostData->hostSparsePermutations;

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
  bulk->volatile_fast_shared_comm_map(stk::topology::NODE_RANK, 0);
  auto& hostVolatileFastSharedCommMapOffset = deviceMeshHostData->hostVolatileFastSharedCommMapOffset;
  auto& hostVolatileFastSharedCommMap = deviceMeshHostData->hostVolatileFastSharedCommMap;

  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::ELEM_RANK; ++rank)
  {
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank].extent(0));
    Kokkos::deep_copy(volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank]);
    Kokkos::deep_copy(volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank]);
  }
}

}
}
