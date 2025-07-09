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

#ifndef STK_MESH_DEVICE_BUCKET_REPOSITORY_HPP
#define STK_MESH_DEVICE_BUCKET_REPOSITORY_HPP

#include "stk_util/stk_config.h"
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/base/DeviceBucket.hpp"
#include "stk_mesh/baseImpl/DevicePartition.hpp"
#include "stk_mesh/baseImpl/NgpMeshImpl.hpp"
#include "Kokkos_Core.hpp"

namespace stk {
namespace mesh {
  
template <typename NgpMemSpace>
class DeviceMeshT;

namespace impl {

template <typename NgpMemSpace>
class DeviceBucketRepository
{
 public:
  using DeviceBucket = DeviceBucketT<NgpMemSpace>;
  using DeviceBucketView = Kokkos::View<DeviceBucket*, stk::ngp::UVMMemSpace>;
  using DeviceBucketUView = Kokkos::View<DeviceBucket*, stk::ngp::UVMMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using DevicePartitionView = Kokkos::View<DevicePartition<NgpMemSpace>*, stk::ngp::UVMMemSpace>;
  using DevicePartitionUView = Kokkos::View<DevicePartition<NgpMemSpace>*, stk::ngp::UVMMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using BoolViewType  = Kokkos::View<bool*, stk::ngp::MemSpace>;

  KOKKOS_DEFAULTED_FUNCTION
  DeviceBucketRepository() = default;

  DeviceBucketRepository(DeviceMeshT<NgpMemSpace>* deviceMesh, [[maybe_unused]] unsigned entityRankCount,
                         unsigned initialBucketCapacity = get_default_initial_bucket_capacity(),
                         unsigned maximumBucketCapacity = get_default_maximum_bucket_capacity())
    : m_mesh(deviceMesh),
      m_initialBucketCapacity(initialBucketCapacity),
      m_maximumBucketCapacity(maximumBucketCapacity)
  {}

  KOKKOS_DEFAULTED_FUNCTION DeviceBucketRepository(DeviceBucketRepository const&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceBucketRepository(DeviceBucketRepository&&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceBucketRepository& operator=(DeviceBucketRepository&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceBucketRepository& operator=(DeviceBucketRepository const&) = default;

  DevicePartitionView& get_partitions(const EntityRank rank)
  {
    return m_partitions[rank];
  }

  DevicePartition<NgpMemSpace>* get_or_create_partition(const EntityRank rank,
                                                        PartOrdinalViewType<NgpMemSpace> const& partOrdinals)
  {
    auto partition = get_partition(rank, partOrdinals);
    if (partition == nullptr) {
      return create_partition(rank, partOrdinals);
    } else {
      return partition;
    }
  }

  DevicePartition<NgpMemSpace>* get_partition(const EntityRank rank, unsigned partitionId)
  {
    STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank.");
    STK_ThrowRequireMsg(partitionId < m_partitions[rank].extent(0), "Invalid partition id.");

    return &(m_partitions[rank](partitionId));
  }
      
  DevicePartition<NgpMemSpace>* get_partition(const EntityRank rank,
                                              PartOrdinalViewType<NgpMemSpace> const& partOrdinals)
  {
    STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank.");

    auto policy = Kokkos::RangePolicy(0, num_partitions(rank));
    DevicePartitionUView compactPartitionUView(m_partitions[rank].data(), m_numPartitions[rank]);
    auto partitionId = search_matching_device_partitions(policy, compactPartitionUView, partOrdinals);

    return (partitionId != INVALID_BUCKET_ID) ? &(m_partitions[rank](partitionId)) : nullptr;
  }

  DevicePartition<NgpMemSpace>* get_partition(DeviceBucket& bucket)
  {
    auto partitionId = bucket.partition_id();
    auto entityRank = bucket.entity_rank();

    return get_partition(entityRank, partitionId);
  }

  DevicePartition<NgpMemSpace>* create_partition(const EntityRank rank,
                                                 PartOrdinalViewType<NgpMemSpace> const& partOrdinals)
  {
    init_or_resize_partition_view(rank);

    DevicePartitionView& resizedPartitions = m_partitions[rank];
    auto endIndex = m_partitionViewEndIdx[rank];

    DevicePartition<NgpMemSpace>* newPartition = new (&resizedPartitions(endIndex)) DevicePartition<NgpMemSpace>(m_mesh, this, partOrdinals, rank, endIndex);
    m_numPartitions[rank]++;
    m_partitionViewEndIdx[rank]++;

    return newPartition;
  }

  DeviceBucket* get_bucket(const EntityRank rank, int bucketId) const
  {
    STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank.");
    STK_ThrowRequireMsg(static_cast<unsigned>(bucketId) < m_buckets[rank].extent(0), "Invalid bucket id.");

    return &(m_buckets[rank](bucketId));
  }

  DeviceBucket* allocate_bucket(const EntityRank rank,
                                PartOrdinalViewType<NgpMemSpace> const& /*partOrdinals*/)
  {
    init_or_resize_bucket_view(rank);

    DeviceBucketView& buckets = m_buckets[rank];
    auto endIndex = m_bucketViewEndIdx[rank];
    DeviceBucket* newBucket = new (&buckets[endIndex]) DeviceBucket();
    newBucket->m_bucketSize = 0;
    newBucket->m_bucketCapacity = get_bucket_capacity();
    newBucket->m_bucketId = endIndex;
    newBucket->m_entityRank = rank;

    m_numBuckets[rank]++;
    m_bucketViewEndIdx[rank]++;
    m_needSyncFromPartitions[rank] = true;

    return newBucket;
  }

  void deallocate_bucket(DeviceBucket* bucket)
  {
    STK_ThrowRequireMsg(bucket != nullptr, "DeviceBucketRepository::deallocate_bucket(): Invalid bucket.");

    auto bucketId = bucket->bucket_id();
    auto bucketRank = bucket->entity_rank();

    m_buckets[bucketRank](bucketId).m_bucketId = INVALID_BUCKET_ID;
    m_needSyncFromPartitions[bucketRank] = true;
    m_numBuckets[bucketRank]--;
  }

  void sync_from_partitions()
  {
    for (auto rank = stk::topology::NODE_RANK; rank < stk::topology::NUM_RANKS; ++rank) {
      sync_from_partitions(rank);
    }
  }

  void sync_from_partitions(EntityRank rank)
  {
    if (m_needSyncFromPartitions[rank]) {
      reset_partition_bucket_count(rank);

      sync_and_sort_bucket_ids(rank);

      reassign_bucket_ids(rank);

      reassign_buckets_owning_partition_ids(rank);

      sync_and_sort_partition_ids(rank);

      sort_buckets_in_partitions(rank);
    }

    // TODO: refactor DeviceMesh first
    // if (m_mesh->get_bulk_on_host().should_sort_buckets_by_first_entity_identifier()) {
    // }

    // TODO: refactor DeviceMesh first
    // update_mesh_indices();
  }

  void reset_partition_bucket_count(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto& partitions = m_partitions[rank]; 
    auto numPartitions = m_numPartitions[rank];

    for (unsigned i = 0; i < numPartitions; ++i) {
      auto& partition = partitions(i);
      partition.reset_bucket_count();

      if (partition.partition_id() == INVALID_PARTITION_ID) {
        --m_numPartitions[rank];
      }
    }
  }

  void sync_and_sort_bucket_ids(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    DeviceBucketUView compactBucketUView(buckets.data(), m_bucketViewEndIdx[rank]);

    if (!Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBucketUView)) {
      Kokkos::sort(compactBucketUView);
    }
    
    // TODO: call removed buckets' destructors here (need to refactor DeviceMesh first)

    m_buckets[rank] = buckets;
    m_bucketViewEndIdx[rank] = m_numBuckets[rank];
  }

  void reassign_bucket_ids(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    for (unsigned i = 0; i < m_numBuckets[rank]; ++i) {
      buckets(i).m_bucketId = i;
    }
  }

  void reassign_buckets_owning_partition_ids(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto numBuckets = m_numBuckets[rank];
    for (unsigned i = 0; i < numBuckets; ++i) {
      auto bucket = buckets(i);
      auto& partition = m_partitions[rank](bucket.m_owningPartitionId);
      partition.add_bucket(&bucket);
    }
  }

  void sync_and_sort_partition_ids(EntityRank rank)
  {
    auto& partitions = m_partitions[rank];
    DevicePartitionUView compactPartitionUView(partitions.data(), m_partitionViewEndIdx[rank]);

    if (!Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactPartitionUView)) {
      Kokkos::sort(compactPartitionUView);
    }
    
    auto endIdx = m_partitionViewEndIdx[rank];
    for (unsigned i = m_numPartitions[rank]; i < endIdx; ++i) {
      if (m_partitions[rank].use_count() == 1) {
        m_partitions[rank](i).~DevicePartition();
      }
    }

    m_partitions[rank] = partitions;
    m_partitionViewEndIdx[rank] = m_numPartitions[rank];
  }

  void sort_buckets_in_partitions(EntityRank rank)
  {
    auto& partitions = m_partitions[rank];
    auto numPartitions = m_numPartitions[rank];
    for (unsigned i = 0; i < numPartitions; ++i) {
      auto& partition = partitions(i);
      partition.reset_partition_id_in_owned_buckets();
      partition.sort_buckets();
    }
  }

  void init_or_resize_partition_view(unsigned rank) { 
    if (m_partitions[rank].extent(0) == 0) {
      m_partitions[rank] = DevicePartitionView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DevicePartitionView"), initialPartitionViewCapacity);
    } else if (m_numPartitions[rank] >= m_partitions[rank].extent(0)) {
      Kokkos::resize(m_partitions[rank], m_partitions[rank].extent(0)*2);
    }
  }

  void init_or_resize_bucket_view(unsigned rank)
  {
    if (m_buckets[rank].extent(0) == 0) {
      m_buckets[rank] = DeviceBucketView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketView"), initialBucketViewCapacity);
    } else if (m_numBuckets[rank] >= m_buckets[rank].extent(0)) {
      Kokkos::resize(m_buckets[rank], m_buckets[rank].extent(0)*2);
    }
  }

  void resize_compact_partition_view(unsigned rank) { 
    Kokkos::resize(m_partitions[rank], m_numPartitions[rank]);
  }

  void resize_compact_bucket_view(unsigned rank) { 
    Kokkos::resize(m_buckets[rank], m_numBuckets[rank]);
  }

  // TODO
  // void update_mesh_indicies() {
  //   FIXME
  //   auto newEntityKeyCount = 0;

  //   Kokkos::resize(Kokkos::WithoutInitializing, m_mesh->hostMeshIndices, newEntityKeyCount);
  //   auto hostMeshIndices = m_mesh->hostMeshIndices;

  //   for (EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::NUM_RANKS; ++rank) {
  //     auto buckets = m_buckets[rank];

  //     for (auto i = 0; i < buckets.extent(0); ++i) {
  //       auto& bucket = buckets(i);
  //       auto id = bucket.bucket_id();

  //       for (unsigned j = 0; j < bucket.size(); ++j) {
  //         hostMeshIndices[ /*offset in bucket */ = FastMeshIndex{id, j};
  //       }
  //     }
  //   }

  //   Kokkos::deep_copy(...);
  // }

  KOKKOS_FUNCTION
  unsigned get_bucket_capacity() const { return m_maximumBucketCapacity; }
  
  KOKKOS_FUNCTION
  unsigned get_initial_bucket_capacity() const { return m_initialBucketCapacity; }

  KOKKOS_FUNCTION
  unsigned get_maximum_bucket_capacity() const { return m_maximumBucketCapacity; }

  size_t num_buckets(EntityRank rank) const { return m_numBuckets[rank]; }
  size_t num_partitions(EntityRank rank) const { return m_numPartitions[rank]; }

  DeviceMeshT<NgpMemSpace>* m_mesh;
  unsigned m_initialBucketCapacity;
  unsigned m_maximumBucketCapacity;

  unsigned m_numBuckets[stk::topology::NUM_RANKS] = {0};
  unsigned m_numPartitions[stk::topology::NUM_RANKS] = {0};

  unsigned m_bucketViewEndIdx[stk::topology::NUM_RANKS] = {0};
  unsigned m_partitionViewEndIdx[stk::topology::NUM_RANKS] = {0};

  DeviceBucketView m_buckets[stk::topology::NUM_RANKS];
  DevicePartitionView m_partitions[stk::topology::NUM_RANKS];
  bool m_needSyncFromPartitions[stk::topology::NUM_RANKS] = {false};
};

} } }


#endif
