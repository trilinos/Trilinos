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

#include <limits>
#include "Kokkos_Macros.hpp"
#include "View/Kokkos_ViewCtor.hpp"
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
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace mesh {
  
template <typename NgpMemSpace>
class DeviceMeshT;

namespace impl {

template <typename NgpMemSpace>
struct DeviceBucketRepoMetadata
{
  using DevPartition = DevicePartition<NgpMemSpace>;
  using DevBucket = DeviceBucketT<NgpMemSpace>;

  static constexpr unsigned INVALID_INDEX = std::numeric_limits<unsigned>::max();

  KOKKOS_FUNCTION
  DeviceBucketRepoMetadata()
  {
    m_deviceNumActiveBuckets = UnsignedViewType<NgpMemSpace>("NumBuckets", stk::topology::NUM_RANKS);

    for (unsigned i = 0; i < stk::topology::END_RANK; ++i) {
      m_bucketViewLastActiveBucketIdx[i] = INVALID_INDEX;
      m_partitionViewLastActivePartitionIdx[i] = INVALID_INDEX;
      m_bucketViewLastInitBucketIdx[i] = INVALID_INDEX;
    }
  }

  void update_repo_meta_bucket_added(DevBucket const& newBucket)
  {
    EntityRank rank = newBucket.entity_rank();
    m_numActiveBuckets[rank]++;
    m_bucketViewLastActiveBucketIdx[rank] = newBucket.bucket_id();
    m_bucketViewLastInitBucketIdx[rank] = newBucket.bucket_id();
  }

  void update_repo_meta_partition_added(DevPartition const& newPartition)
  {
    EntityRank rank = newPartition.get_rank();
    m_numActivePartitions[rank]++;
    m_partitionViewLastActivePartitionIdx[rank] = newPartition.partition_id();
  }

  void update_repo_meta_bucket_removed(DevBucket const& removedBucket)
  {
    EntityRank rank = removedBucket.entity_rank();
    m_numActiveBuckets[rank]--;
  }

  void update_repo_meta_partition_removed(DevPartition const& removedPartition)
  {
    EntityRank rank = removedPartition.get_rank();
    m_numActivePartitions[rank]--;
  }

  KOKKOS_FUNCTION
  unsigned num_buckets(EntityRank rank) const {
    KOKKOS_IF_ON_HOST((
      return m_numActiveBuckets[rank];
    ))
    KOKKOS_IF_ON_DEVICE((
      return m_deviceNumActiveBuckets(rank);
    ))
  }

  KOKKOS_FUNCTION
  unsigned num_partitions(EntityRank rank) const {
    KOKKOS_IF_ON_HOST((
      return m_numActivePartitions[rank];
    ))
    KOKKOS_IF_ON_DEVICE((
      return m_deviceNumActivePartitions(rank);
    ))
  }

  unsigned get_last_active_bucket_idx(EntityRank rank) const {
    return m_bucketViewLastActiveBucketIdx[rank];
  }

  unsigned get_last_active_partition_idx(EntityRank rank) const {
    return m_partitionViewLastActivePartitionIdx[rank];
  }

  unsigned get_active_bucket_span(EntityRank rank) const
  {
    auto lastBucketIdx = m_bucketViewLastActiveBucketIdx[rank];
    auto modBucketCount = (lastBucketIdx == INVALID_INDEX) ? 0 : lastBucketIdx+1;
    return modBucketCount;
  }

  unsigned get_active_partition_span(EntityRank rank) const
  {
    auto lastPartitionIdx = m_partitionViewLastActivePartitionIdx[rank];
    auto modPartitionCount = (lastPartitionIdx == INVALID_INDEX) ? 0 : lastPartitionIdx+1;
    return modPartitionCount;
  }

  unsigned get_next_avail_bucket_idx(EntityRank rank) const
  {
    auto lastBucketIdx = m_bucketViewLastActiveBucketIdx[rank];
    auto nextBucketIdx = (lastBucketIdx == INVALID_INDEX) ? 0 : lastBucketIdx+1;
    return nextBucketIdx;
  }

  unsigned get_next_avail_partition_idx(EntityRank rank) const
  {
    auto lastPartitionIdx = m_partitionViewLastActivePartitionIdx[rank];
    auto nextPartitionIdx = (lastPartitionIdx == INVALID_INDEX) ? 0 : lastPartitionIdx+1;
    return nextPartitionIdx;
  }

  void set_to_match_active_buckets(EntityRank rank)
  {
    STK_ThrowRequire(m_numActiveBuckets[rank] > 0);
    m_bucketViewLastActiveBucketIdx[rank] = m_numActiveBuckets[rank]-1;
  }

  void set_to_match_active_partitions(EntityRank rank)
  {
    STK_ThrowRequire(m_numActivePartitions[rank] > 0);
    m_partitionViewLastActivePartitionIdx[rank] = m_numActivePartitions[rank]-1;
  }

  bool need_sync_from_partition(EntityRank rank) const { return m_needSyncFromPartitions[rank]; }

  void set_need_sync_from_partition(EntityRank rank, bool needSync) {
    m_needSyncFromPartitions[rank] = m_needSyncFromPartitions[rank] | needSync;
  }

  void sync_num_buckets_host_to_device()
  {
    auto hostNumActiveBuckets = Kokkos::create_mirror_view(m_deviceNumActiveBuckets);
    for (unsigned i = 0; i < hostNumActiveBuckets.extent(0); ++i) {
      hostNumActiveBuckets(i) = m_numActiveBuckets[i];
    }
    Kokkos::deep_copy(m_deviceNumActiveBuckets, hostNumActiveBuckets);
  }

  void sync_num_partitions_host_to_device()
  {
    auto hostNumActivePartitions = Kokkos::create_mirror_view(m_deviceNumActivePartitions);
    for (unsigned i = 0; i < hostNumActivePartitions.extent(0); ++i) {
      hostNumActivePartitions(i) = m_numActivePartitions[i];
    }
    Kokkos::deep_copy(m_deviceNumActivePartitions, hostNumActivePartitions);
  }

  void reset_bucket_counts(EntityRank rank)
  {
    m_numActiveBuckets[rank] = 0;
    m_bucketViewLastActiveBucketIdx[rank] = INVALID_INDEX;
  }

  void reset_partition_counts(EntityRank rank)
  {
    m_numActivePartitions[rank] = 0;
    m_partitionViewLastActivePartitionIdx[rank] = INVALID_INDEX;
  }

  UnsignedViewType<NgpMemSpace> m_deviceNumActiveBuckets;
  UnsignedViewType<NgpMemSpace> m_deviceNumActivePartitions;
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_numActiveBuckets{0};
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_numActivePartitions{0};
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_bucketViewLastActiveBucketIdx{0};
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_partitionViewLastActivePartitionIdx{0};
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_bucketViewLastInitBucketIdx{0};
  Kokkos::Array<bool, stk::topology::NUM_RANKS> m_needSyncFromPartitions{false};
};

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

  DeviceBucketRepository(DeviceMeshT<NgpMemSpace>* deviceMesh,
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

  KOKKOS_FUNCTION
  ~DeviceBucketRepository() {
    clear_device_buckets_and_partitions();
  }

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
    Kokkos::Profiling::pushRegion("get_partition");
    STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank.");

    auto policy = Kokkos::RangePolicy(0, get_active_partition_span(rank));
    DevicePartitionUView compactPartitionUView(m_partitions[rank].data(), get_active_partition_span(rank));
    auto partitionId = search_matching_device_partitions(policy, compactPartitionUView, partOrdinals);
    Kokkos::Profiling::popRegion();

    return (partitionId != INVALID_BUCKET_ID) ? &(m_partitions[rank](partitionId)) : nullptr;
  }

  DevicePartition<NgpMemSpace>* get_partition(DeviceBucket& bucket)
  {
    auto partitionId = bucket.partition_id();
    auto entityRank = bucket.entity_rank();

    return get_partition(entityRank, partitionId);
  }

  template <typename PartOrdinals>
  DevicePartition<NgpMemSpace>* create_partition(const EntityRank rank, PartOrdinals const& partOrdinals)
  {
    Kokkos::Profiling::pushRegion("create_partition");
    init_or_expand_partition_view(rank);

    DevicePartitionView& resizedPartitions = m_partitions[rank];
    auto newPartitionIdx = get_next_avail_partition_idx(rank);

    new (&resizedPartitions(newPartitionIdx)) DevicePartition<NgpMemSpace>(m_mesh, this, partOrdinals, rank, newPartitionIdx);
    auto& newPartition = resizedPartitions(newPartitionIdx);

    m_repoMeta.update_repo_meta_partition_added(newPartition);
    Kokkos::Profiling::popRegion();

    return &newPartition;
  }

  DeviceBucket* get_bucket(const EntityRank rank, int bucketId) const
  {
    STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank.");
    STK_ThrowRequireMsg(static_cast<unsigned>(bucketId) < m_buckets[rank].extent(0), "Invalid bucket id.");

    return &(m_buckets[rank](bucketId));
  }

  DeviceBucket* construct_new_bucket_without_part_ordinals(const EntityRank rank)
  {
    init_or_expand_bucket_view(rank);

    DeviceBucketView& buckets = m_buckets[rank];
    auto newBucketIdx = get_next_avail_bucket_idx(rank);

    new (&buckets[newBucketIdx]) DeviceBucket();
    DeviceBucket& newBucket = buckets[newBucketIdx];

    newBucket.m_bucketSize = 0;
    newBucket.m_bucketCapacity = get_bucket_capacity();
    newBucket.m_bucketId = newBucketIdx;
    newBucket.m_entityRank = rank;

    m_repoMeta.update_repo_meta_bucket_added(newBucket);
    m_repoMeta.set_need_sync_from_partition(rank, true);
    
    return &newBucket;
  }

  DeviceBucket* construct_new_bucket(const EntityRank rank,
                                     PartOrdinalViewType<NgpMemSpace> const& partOrdinals)
  {
    Kokkos::Profiling::pushRegion("construct_bucket");
    auto newBucket = construct_new_bucket_without_part_ordinals(rank);
    newBucket->m_partOrdinals = partOrdinals;
    Kokkos::Profiling::popRegion();
    return newBucket;
  }

  void invalidate_bucket(DeviceBucket* bucket)
  {
    STK_ThrowRequireMsg(bucket != nullptr, "DeviceBucketRepository::invalidate_bucket(): Invalid bucket pointer.");

    auto bucketId = bucket->bucket_id();
    auto bucketRank = bucket->entity_rank();

    auto& removedBucket = m_buckets[bucketRank](bucketId);
    removedBucket.m_bucketId = INVALID_BUCKET_ID;
    m_repoMeta.update_repo_meta_bucket_removed(removedBucket);
    m_repoMeta.set_need_sync_from_partition(bucketRank, true);
  }

  void sync_num_buckets_host_to_device() {
    m_repoMeta.sync_num_buckets_host_to_device();
  }

  void sync_num_partitions_host_to_device() {
    m_repoMeta.sync_num_partitions_host_to_device();
  }

  void sync_from_partitions()
  {
    for (auto rank = stk::topology::NODE_RANK; rank < stk::topology::NUM_RANKS; ++rank) {
      sync_from_partitions(rank);
    }
    sync_num_buckets_host_to_device();
    sync_num_partitions_host_to_device();
  }

  void sync_from_partitions(EntityRank rank)
  {
    Kokkos::Profiling::pushRegion("sync_from_partiitons");
    if (m_repoMeta.need_sync_from_partition(rank)) {
      reset_partition_bucket_count(rank);

      sync_and_sort_bucket_ids(rank);

      reassign_bucket_ids(rank);

      reassign_buckets_owning_partition_ids(rank);

      sync_and_sort_partition_ids(rank);

      sort_buckets_in_partitions(rank);
    }
    Kokkos::Profiling::popRegion();

    // TODO
    // if (m_mesh->get_bulk_on_host().should_sort_buckets_by_first_entity_identifier()) {
    // }

    // TODO
    // update_buckets_mesh_indices(rank);
  }

  void sync_and_sort_bucket_ids(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    DeviceBucketUView compactBucketUView(buckets.data(), get_active_bucket_span(rank));

    if (!Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactBucketUView)) {
      Kokkos::sort(compactBucketUView);
    }
    
    m_buckets[rank] = buckets;
    m_repoMeta.set_to_match_active_buckets(rank);
  }

  void reassign_bucket_ids(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto numActiveBuckets = num_buckets(rank);
    for (unsigned i = 0; i < numActiveBuckets; ++i) {
      buckets(i).m_bucketId = i;
    }
  }

  void reassign_buckets_owning_partition_ids(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto numActiveBuckets = num_buckets(rank);
    for (unsigned i = 0; i < numActiveBuckets; ++i) {
      auto bucket = buckets(i);
      auto& partition = m_partitions[rank](bucket.m_owningPartitionId);
      partition.add_bucket(&bucket);
    }
  }

  void sync_and_sort_partition_ids(EntityRank rank)
  {
    auto& partitions = m_partitions[rank];
    DevicePartitionUView compactPartitionUView(partitions.data(), get_active_partition_span(rank));

    if (!Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactPartitionUView)) {
      Kokkos::sort(compactPartitionUView);
    }
    
    m_partitions[rank] = partitions;
    m_repoMeta.set_to_match_active_partitions(rank);
  }

  void sort_buckets_in_partitions(EntityRank rank)
  {
    auto& partitions = m_partitions[rank];
    auto numPartitions = num_partitions(rank);
    for (unsigned i = 0; i < numPartitions; ++i) {
      auto& partition = partitions(i);
      partition.reset_partition_id_in_owned_buckets();
      partition.sort_buckets();
    }
  }

  void init_or_expand_partition_view(EntityRank rank, unsigned hostPartitionVectorCapcity = initialDevicePartitionViewCapacity)
  { 
    Kokkos::Profiling::pushRegion("init_or_expand_partition_view");
    if (m_partitions[rank].extent(0) == 0) {
      m_partitions[rank] = DevicePartitionView(Kokkos::view_alloc("DevicePartitionView", Kokkos::WithoutInitializing), hostPartitionVectorCapcity);
    } else if (get_next_avail_partition_idx(rank) >= m_partitions[rank].extent(0)) {
      Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_partitions[rank], m_partitions[rank].extent(0)*2);
    }
    Kokkos::Profiling::popRegion();
  }

  void init_or_expand_bucket_view(EntityRank rank, unsigned hostBucketVectorCapacity = initialDeviceBucketViewCapacity)
  {
    Kokkos::Profiling::pushRegion("init_or_expand_bucket_view");
    if (m_buckets[rank].extent(0) == 0) {
      m_buckets[rank] = DeviceBucketView(Kokkos::view_alloc("DeviceBucketView", Kokkos::WithoutInitializing), hostBucketVectorCapacity);
    } else if (get_next_avail_bucket_idx(rank) >= m_buckets[rank].extent(0)) {
      Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_buckets[rank], m_buckets[rank].extent(0)*2);
    }
    Kokkos::Profiling::popRegion();
  }

  // void resize_compact_partition_view(EntityRank rank) { 
  //   // Kokkos::resize(Kokkos::view_alloc(Kokkos::SequentialHostInit), m_partitions[rank], m_numActivePartitions[rank]);
  //   Kokkos::resize(m_partitions[rank], m_numActivePartitions[rank]);
  //   // FIXME destroy removed partitions
  // }

  // void resize_compact_bucket_view(EntityRank rank) { 
  //   // Kokkos::resize(Kokkos::view_alloc(Kokkos::SequentialHostInit), m_buckets[rank], m_numActiveBuckets[rank]);
  //   Kokkos::resize(m_buckets[rank], m_numActiveBuckets[rank]);
  //   // FIXME destroy removed buckets
  // }

  void reset_partition_bucket_count(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto& partitions = m_partitions[rank]; 
    auto numPartitions = num_partitions(rank);

    for (unsigned i = 0; i < numPartitions; ++i) {
      auto& partition = partitions(i);
      partition.reset_bucket_count();

      if (partition.partition_id() == INVALID_PARTITION_ID) {
        m_repoMeta.update_repo_meta_partition_removed(partition);
      }
    }
  }

  void reset_and_invalidate_all_buckets(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto bucketCount = get_active_bucket_span(rank);

    for (unsigned i = 0; i < bucketCount; ++i) {
      auto& bucket = buckets(i);
      bucket.m_bucketId = INVALID_BUCKET_ID;
    }
    m_repoMeta.reset_bucket_counts(rank);
  }

  void reset_and_invalidate_all_partitions(EntityRank rank)
  {
    auto& partitions = m_partitions[rank];
    auto partitionCount = get_active_partition_span(rank);

    for (unsigned i = 0; i < partitionCount; ++i) {
      auto& partition = partitions(i);
      partition.m_partitionId = INVALID_PARTITION_ID;
    }
    m_repoMeta.reset_partition_counts(rank);
  }

  void reset_and_invalidate_all(EntityRank rank)
  {
    reset_and_invalidate_all_buckets(rank);
    reset_and_invalidate_all_partitions(rank);
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

  KOKKOS_FUNCTION
  unsigned num_buckets(EntityRank rank) const {
    return m_repoMeta.num_buckets(rank);
  }

  KOKKOS_FUNCTION
  unsigned num_partitions(EntityRank rank) const {
    return m_repoMeta.num_partitions(rank);
  }

  unsigned get_last_active_bucket_idx(EntityRank rank) const {
    return m_repoMeta.get_last_active_bucket_idx(rank);
  }

  unsigned get_last_active_partition_idx(EntityRank rank) const {
    return m_repoMeta.get_last_active_partition_idx(rank);
  }

  unsigned get_active_bucket_span(EntityRank rank) const {
    return m_repoMeta.get_active_bucket_span(rank);
  }

  unsigned get_active_partition_span(EntityRank rank) const {
    return m_repoMeta.get_active_partition_span(rank);
  }

  unsigned get_next_avail_bucket_idx(EntityRank rank) const
  {
    return m_repoMeta.get_next_avail_bucket_idx(rank);
  }

  unsigned get_next_avail_partition_idx(EntityRank rank) const {
    return m_repoMeta.get_next_avail_partition_idx(rank);
  }

  KOKKOS_FUNCTION
  bool is_last_partition_reference(unsigned rank = stk::topology::NODE_RANK) const {
    return (m_partitions[rank].use_count() == 1);
  }

  KOKKOS_FUNCTION
  bool is_last_bucket_reference(unsigned rank = stk::topology::NODE_RANK) const {
    return (m_buckets[rank].use_count() == 1);
  }

  KOKKOS_FUNCTION
  void clear_device_buckets_and_partitions()
  {
    KOKKOS_IF_ON_HOST((
      for (EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::END_RANK; rank++) {
        if (is_last_partition_reference(rank)) {
          auto partitions = std::max(get_active_partition_span(rank), num_partitions(rank));
          for (unsigned iPartition = 0; iPartition < partitions; ++iPartition) {
            m_partitions[rank](iPartition).~DevicePartition();
          }
          if (is_last_bucket_reference(rank)) {
            auto buckets = std::max(get_active_bucket_span(rank), num_buckets(rank));
            for (unsigned iBucket = 0; iBucket < buckets; ++iBucket) {
              m_buckets[rank](iBucket).~DeviceBucket();
            }
          }
        }
      }
    ))
  }

  // copy from host
  template <typename HostBuckets, typename HostPartitions>
  void copy_buckets_and_partitions_from_host(EntityRank rank, HostBuckets const& hostBuckets, HostPartitions const& hostPartitions, bool& anyBucketChanges)
  {
    Kokkos::Profiling::pushRegion("copy_buckets_and_partitions_from_host");
    auto numHostPartitions = hostPartitions.size();

    if (hostBuckets.size() != num_buckets(rank)) {
      anyBucketChanges = true;
    }

    reset_and_invalidate_all(rank);

    auto bucketCapacityToInit = hostBuckets.size();
    auto partitionCapacityToInit = hostPartitions.size();
    DeviceBucketView bucketBuffer(Kokkos::view_alloc(Kokkos::WithoutInitializing, "BucketBuffer"), bucketCapacityToInit);
    DevicePartitionView partitionBuffer(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PartitionBuffer"), partitionCapacityToInit);

    for (unsigned partitionIdx = 0, bucketIdxInView = 0; partitionIdx < numHostPartitions; ++partitionIdx) {
      auto hostPartition = hostPartitions[partitionIdx];
      auto& partitionParts = hostPartition->get_legacy_partition_id();

      new (&partitionBuffer[partitionIdx]) impl::DevicePartition<NgpMemSpace>(m_mesh, this, partitionParts, rank, partitionIdx);

      // Init partition inner bucketPtr view
      partitionBuffer[partitionIdx].init_bucket_view(hostPartition->num_buckets());

      auto numBucketsInPartition = hostPartition->num_buckets();
      for (unsigned i = 0; i < numBucketsInPartition; ++i, ++bucketIdxInView) {
        auto& iBucket = bucketIdxInView;
        stk::mesh::Bucket& stkBucket = *hostBuckets[iBucket];
        const int ngpBucketId = stkBucket.ngp_mesh_bucket_id();

        if (ngpBucketId == INVALID_BUCKET_ID) {
          new (&bucketBuffer[iBucket]) DeviceBucketT<NgpMemSpace>();
          bucketBuffer[iBucket].initialize_bucket_attributes(stkBucket);
          bucketBuffer[iBucket].initialize_part_ordinals_from_host(stkBucket);
          bucketBuffer[iBucket].update_entity_data_from_host(stkBucket);
          bucketBuffer[iBucket].update_sparse_connectivity_from_host(stkBucket);
          anyBucketChanges = true;
        } else {
          new (&bucketBuffer[iBucket]) DeviceBucketT(m_buckets[rank][ngpBucketId]);
          if (stkBucket.ngp_mesh_bucket_is_modified()) {
            bucketBuffer[iBucket].update_entity_data_from_host(stkBucket);
            bucketBuffer[iBucket].update_sparse_connectivity_from_host(stkBucket);
            anyBucketChanges = true;
          }
          bucketBuffer[iBucket].m_bucketId = stkBucket.bucket_id();
        }

        stkBucket.set_ngp_mesh_bucket_id(iBucket);

        // Update Device Bucket View metadata
        STK_ThrowRequire(iBucket == stkBucket.bucket_id());
        m_repoMeta.update_repo_meta_bucket_added(bucketBuffer[iBucket]);

        // Update device buckets in partition
        STK_ThrowRequire(i == partitionBuffer[partitionIdx].get_next_avail_bucket_idx());
        partitionBuffer[partitionIdx].add_bucket(&bucketBuffer[iBucket]);
      }

      // Update Device Partition View metadata
      m_repoMeta.update_repo_meta_partition_added(partitionBuffer[partitionIdx]);
    }

    if (is_last_bucket_reference(rank)) {
      for (unsigned i = 0; i < m_buckets[rank].size(); ++i) {
        m_buckets[rank][i].~DeviceBucketT();
      }
    }

    if (is_last_partition_reference(rank)) {
      for (unsigned i = 0; i < m_partitions[rank].size(); ++i) {
        m_partitions[rank][i].~DevicePartition();
      }
    }

    m_buckets[rank] = bucketBuffer;
    m_partitions[rank] = partitionBuffer;
    Kokkos::Profiling::popRegion();
  }

  unsigned m_initialBucketCapacity;
  unsigned m_maximumBucketCapacity;

  DeviceMeshT<NgpMemSpace>* m_mesh;
  DeviceBucketRepoMetadata<NgpMemSpace> m_repoMeta;
  
  DeviceBucketView m_buckets[stk::topology::NUM_RANKS];
  DevicePartitionView m_partitions[stk::topology::NUM_RANKS];
};

} } }


#endif
