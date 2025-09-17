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
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Macros.hpp"
#include "View/Kokkos_ViewCtor.hpp"
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
  using BoolViewType = Kokkos::View<bool*, ngp::UVMMemSpace>;

  static constexpr unsigned INVALID_INDEX = std::numeric_limits<unsigned>::max();

  DeviceBucketRepoMetadata()
  {
    m_deviceNumActiveBuckets = UnsignedViewType<NgpMemSpace>("NumBuckets", stk::topology::NUM_RANKS);
    m_deviceNumActivePartitions = UnsignedViewType<NgpMemSpace>("NumPartitions", stk::topology::NUM_RANKS);
    m_needSyncFromPartitions = BoolViewType("NeedSyncFromDevice", stk::topology::NUM_RANKS);

    for (unsigned i = 0; i < stk::topology::END_RANK; ++i) {
      m_bucketViewLastActiveBucketIdx[i] = INVALID_INDEX;
      m_partitionViewLastActivePartitionIdx[i] = INVALID_INDEX;
      m_bucketViewLastInitBucketIdx[i] = INVALID_INDEX;
    }
  }

  KOKKOS_DEFAULTED_FUNCTION
  ~DeviceBucketRepoMetadata() = default;

  void update_repo_meta_bucket_added(DevBucket const& newBucket)
  {
    EntityRank rank = newBucket.entity_rank();
    STK_ThrowRequireMsg(rank!=stk::topology::INVALID_RANK,"update_repo_meta_bucket_added: invalid rank="<<rank);
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

  KOKKOS_INLINE_FUNCTION
  unsigned num_buckets(EntityRank rank) const {
    KOKKOS_IF_ON_HOST((
      return m_numActiveBuckets[rank];
    ))
    KOKKOS_IF_ON_DEVICE((
      return m_deviceNumActiveBuckets(rank);
    ))
  }

  KOKKOS_INLINE_FUNCTION
  unsigned num_partitions(EntityRank rank) const {
    KOKKOS_IF_ON_HOST((
      return m_numActivePartitions[rank];
    ))
    KOKKOS_IF_ON_DEVICE((
      return m_deviceNumActivePartitions(rank);
    ))
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_last_active_bucket_idx(EntityRank rank) const {
    return m_bucketViewLastActiveBucketIdx[rank];
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_last_active_partition_idx(EntityRank rank) const {
    return m_partitionViewLastActivePartitionIdx[rank];
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_active_bucket_span(EntityRank rank) const
  {
    auto lastBucketIdx = m_bucketViewLastActiveBucketIdx[rank];
    auto modBucketCount = (lastBucketIdx == INVALID_INDEX || num_buckets(rank) == 0) ? 0 : lastBucketIdx+1;
    return modBucketCount;
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_active_partition_span(EntityRank rank) const
  {
    auto lastPartitionIdx = m_partitionViewLastActivePartitionIdx[rank];
    auto modPartitionCount = (lastPartitionIdx == INVALID_INDEX || num_partitions(rank) == 0) ? 0 : lastPartitionIdx+1;
    return modPartitionCount;
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_next_avail_bucket_idx(EntityRank rank) const
  {
    auto lastBucketIdx = m_bucketViewLastActiveBucketIdx[rank];
    auto nextBucketIdx = (lastBucketIdx == INVALID_INDEX) ? 0 : lastBucketIdx+1;
    return nextBucketIdx;
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_next_avail_partition_idx(EntityRank rank) const
  {
    auto lastPartitionIdx = m_partitionViewLastActivePartitionIdx[rank];
    auto nextPartitionIdx = (lastPartitionIdx == INVALID_INDEX) ? 0 : lastPartitionIdx+1;
    return nextPartitionIdx;
  }

  void set_to_match_active_buckets(EntityRank rank)
  {
    m_bucketViewLastActiveBucketIdx[rank] = (m_numActiveBuckets[rank] == 0) ? INVALID_INDEX
                                                                            : m_numActiveBuckets[rank]-1;
  }

  void set_to_match_active_partitions(EntityRank rank)
  {
    m_partitionViewLastActivePartitionIdx[rank] = (m_numActivePartitions[rank] == 0) ? INVALID_INDEX
                                                                                     : m_numActivePartitions[rank]-1;
  }

  KOKKOS_INLINE_FUNCTION
  bool need_sync_from_partition(EntityRank rank) const { return m_needSyncFromPartitions[rank]; }

  KOKKOS_INLINE_FUNCTION
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
  BoolViewType m_needSyncFromPartitions;

  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_numActiveBuckets{0};
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_numActivePartitions{0};
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_bucketViewLastActiveBucketIdx{0};
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_partitionViewLastActivePartitionIdx{0};
  Kokkos::Array<unsigned, stk::topology::NUM_RANKS> m_bucketViewLastInitBucketIdx{0};
};

struct PartOrdinalsProxyIndices {
  KOKKOS_DEFAULTED_FUNCTION
  PartOrdinalsProxyIndices() = default;

  KOKKOS_INLINE_FUNCTION
  PartOrdinalsProxyIndices(EntityRank entityRank, Ordinal* ptr, unsigned len)
    : rank(entityRank),
      startPtr(ptr),
      length(len)
  {}

  KOKKOS_DEFAULTED_FUNCTION
  PartOrdinalsProxyIndices(PartOrdinalsProxyIndices const& other) = default;

  KOKKOS_INLINE_FUNCTION
  bool operator==(PartOrdinalsProxyIndices const& other) const
  {
    if (rank != other.rank) { return false; }
    if (length != other.length) { return false; }

    for (unsigned i = 0; i < length; ++i) {
      auto thisOrd = startPtr + i;
      auto otherOrd = other.startPtr + i;
      if (*thisOrd != *otherOrd) { return false; }
    }
    return true;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator<(PartOrdinalsProxyIndices const& other) const
  {
    if (rank != other.rank) { return rank < other.rank; }
    if (length != other.length) { return length < other.length; }

    for (unsigned i = 0; i < length; ++i) {
      auto thisOrd = startPtr + i;
      auto otherOrd = other.startPtr + i;
      if (*thisOrd != *otherOrd) {
        return *thisOrd < *otherOrd;
      }
    }
    return false;
  }

  EntityRank rank;
  Ordinal* startPtr;
  unsigned length;
};

template <typename NgpMemSpace>
class DeviceBucketRepository
{
 public:
  using DeviceMesh = DeviceMeshT<NgpMemSpace>;
  using DeviceBucket = DeviceBucketT<NgpMemSpace>;
  using DeviceBucketView = Kokkos::View<DeviceBucket*, stk::ngp::UVMMemSpace>;
  using DeviceBucketUView = Kokkos::View<DeviceBucket*, stk::ngp::UVMMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using DeviceBucketWrapperView = Kokkos::View<impl::DeviceBucketWrapper<NgpMemSpace>*, stk::ngp::UVMMemSpace>;
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

  KOKKOS_INLINE_FUNCTION
  ~DeviceBucketRepository() {
    clear_device_buckets_and_partitions();
  }

  KOKKOS_INLINE_FUNCTION
  DevicePartitionView& get_partitions(const EntityRank rank)
  {
    return m_partitions[rank];
  }

  DevicePartition<NgpMemSpace>* get_or_create_partition(const EntityRank rank,
                                                        PartOrdinalViewType<NgpMemSpace> const& partOrdinals,
                                                        bool copyPartOrdinals = false)
  {
    auto partition = get_partition(rank, partOrdinals);
    if (partition == nullptr) {
      return create_partition(rank, partOrdinals, copyPartOrdinals);
    } else {
      return partition;
    }
  }

  KOKKOS_INLINE_FUNCTION
  DevicePartition<NgpMemSpace>* get_partition(const EntityRank rank, unsigned partitionId) const
  {
    KOKKOS_IF_ON_HOST((
      STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank="<<rank);
      STK_ThrowRequireMsg(partitionId < m_partitions[rank].extent(0), "Invalid partition id.");
    ))

    return &(m_partitions[rank](partitionId));
  }
      
  DevicePartition<NgpMemSpace>* get_partition(const EntityRank rank,
                                              PartOrdinalViewType<NgpMemSpace> const& partOrdinals) const
  {
    Kokkos::Profiling::pushRegion("get_partition");
    STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank="<<rank);

    auto policy = Kokkos::RangePolicy(0, get_active_partition_span(rank));
    DevicePartitionUView compactPartitionUView(m_partitions[rank].data(), get_active_partition_span(rank));
    auto partitionId = search_matching_device_partitions(policy, compactPartitionUView, partOrdinals);
    Kokkos::Profiling::popRegion();

    return (partitionId != INVALID_BUCKET_ID) ? &(m_partitions[rank](partitionId)) : nullptr;
  }

  DevicePartition<NgpMemSpace>* get_partition(DeviceBucket& bucket) const
  {
    auto partitionId = bucket.partition_id();
    auto entityRank = bucket.entity_rank();

    return get_partition(entityRank, partitionId);
  }

  template <typename EntityViewType, typename PartOrdinalProxyIndicesViewType>
  void batch_get_partitions(EntityViewType const& entities, PartOrdinalProxyIndicesViewType const& partOrdinalProxy, UnsignedPairViewType<NgpMemSpace>& srcDestPartitionIds) const;

  template <typename PartOrdinals>
  DevicePartition<NgpMemSpace>* create_partition(const EntityRank rank, PartOrdinals const& partOrdinals, bool copyPartOrdinals = false)
  {
    Kokkos::Profiling::pushRegion("create_partition");
    init_or_expand_partition_view(rank);

    DevicePartitionView& resizedPartitions = m_partitions[rank];
    auto newPartitionIdx = get_next_avail_partition_idx(rank);

    new (&resizedPartitions(newPartitionIdx)) DevicePartition<NgpMemSpace>(m_mesh, this, partOrdinals, rank, newPartitionIdx, copyPartOrdinals);
    auto& newPartition = resizedPartitions(newPartitionIdx);

    m_repoMeta.update_repo_meta_partition_added(newPartition);
    Kokkos::Profiling::popRegion();

    return &newPartition;
  }

  template <typename EntityViewType>
  void batch_move_entities(EntityViewType const& entities, UnsignedPairViewType<NgpMemSpace> const& srcDestPartitions)
  {
    UnsignedViewType<NgpMemSpace> srcBucketIds("srcBucketIds", entities.extent(0));
    auto& deviceMesh = *m_mesh;

    Kokkos::parallel_for(entities.extent(0),
      KOKKOS_CLASS_LAMBDA(const int i) {
        auto entity = entities(i);
        auto rank = deviceMesh.entity_rank(entity);
        auto srcBucketId = deviceMesh.fast_mesh_index(entity).bucket_id;
        auto srcPartitionId = srcDestPartitions(i).first;
        auto& partition = m_partitions[rank](srcPartitionId);

        srcBucketIds(i) = srcBucketId;
        partition.remove_entity(deviceMesh, entity);
      }
    );

    Kokkos::fence();

    auto hostEntities = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace{}, entities);
    auto hostSrcDestPartitions = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace{}, srcDestPartitions);
    auto hostSrcBucketIds = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace{}, srcBucketIds);

    for (unsigned i = 0; i < hostEntities.extent(0); ++i) {
      auto entity = hostEntities(i);
      auto rank = m_mesh->get_bulk_on_host().entity_rank(entity);
      auto destPartitionId = hostSrcDestPartitions(i).second;
      auto destPartition = get_partition(rank, destPartitionId);
      auto srcBucketId = hostSrcBucketIds(i);

      destPartition->add_entity(entity, srcBucketId);
      set_need_sync_from_partition(rank, true);
    }

    Kokkos::fence();
  }

  template <typename PartOrdinalIndicesViewType>
  void batch_create_partitions(PartOrdinalIndicesViewType const& aggregatedPartOrdinalsProxy)
  {
    using UView = Kokkos::View<Ordinal*, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    auto partOrdinalsProxy = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), aggregatedPartOrdinalsProxy);

    for (unsigned i = 0; i < partOrdinalsProxy.extent(0); ++i) {
      auto rank = partOrdinalsProxy(i).rank;
      auto startPtr = partOrdinalsProxy(i).startPtr;
      auto length = partOrdinalsProxy(i).length;

      if (startPtr == nullptr) { continue; }

      UView view(startPtr, length);
      get_or_create_partition(rank, view, true);
    }
  }

  KOKKOS_INLINE_FUNCTION
  DeviceBucket* get_bucket(const EntityRank rank, int bucketIdx) const
  {
    KOKKOS_IF_ON_HOST((
      STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank.");
      STK_ThrowRequireMsg(static_cast<unsigned>(bucketIdx) < m_buckets[rank].extent(0), "Invalid bucket id.");
    ))

    return &(m_buckets[rank](bucketIdx));
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
    set_need_sync_from_partition(rank, true);
    
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

  void copy_bucket_connectivity(const EntityRank rank, unsigned srcBucketId, unsigned destBucketId)
  {
    auto srcBucket = get_bucket(rank, srcBucketId);
    auto destBucket = get_bucket(rank, destBucketId);

    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), destBucket->m_nodeConnectivity, srcBucket->m_nodeConnectivity.extent(0));
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), destBucket->m_nodeConnectivityOffsets, srcBucket->m_nodeConnectivityOffsets.extent(0));
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), destBucket->m_nodeOrdinals, srcBucket->m_nodeOrdinals.extent(0));
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), destBucket->m_sparseConnectivityOffsets, srcBucket->m_sparseConnectivityOffsets.extent(0), srcBucket->m_sparseConnectivityOffsets.extent(1));
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), destBucket->m_sparseConnectivity, srcBucket->m_sparseConnectivity.extent(0));
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), destBucket->m_sparseConnectivityOrdinals, srcBucket->m_sparseConnectivityOrdinals.extent(0));

    auto execSpace = typename NgpMemSpace::execution_space{};
    Kokkos::Experimental::copy(execSpace, srcBucket->m_nodeConnectivity, destBucket->m_nodeConnectivity);
    Kokkos::Experimental::copy(execSpace, srcBucket->m_nodeConnectivityOffsets, destBucket->m_nodeConnectivityOffsets);
    Kokkos::Experimental::copy(execSpace, srcBucket->m_nodeOrdinals, destBucket->m_nodeOrdinals);
    Kokkos::deep_copy(destBucket->m_sparseConnectivityOffsets, srcBucket->m_sparseConnectivityOffsets);
    Kokkos::Experimental::copy(execSpace, srcBucket->m_sparseConnectivity, destBucket->m_sparseConnectivity);
    Kokkos::Experimental::copy(execSpace, srcBucket->m_sparseConnectivityOrdinals, destBucket->m_sparseConnectivityOrdinals);
  }

  void invalidate_bucket(DeviceBucket* bucket)
  {
    STK_ThrowRequireMsg(bucket != nullptr, "Invalid bucket pointer.");
    STK_ThrowRequireMsg(bucket->size() == 0, "Bucket is not empty.");

    bucket->m_bucketId = INVALID_BUCKET_ID;
    auto partition = get_partition(bucket->entity_rank(), bucket->partition_id());
    partition->update_partition_meta_bucket_removed(bucket);
    m_repoMeta.update_repo_meta_bucket_removed(*bucket);
    set_need_sync_from_partition(bucket->entity_rank(), true);
  }

  void invalidate_partition(DevicePartition<NgpMemSpace>* partition)
  {
    STK_ThrowRequireMsg(partition != nullptr, "Invalid partition pointer.");
    STK_ThrowRequireMsg(partition->num_buckets() == 0, "Partition is not empty.");
    STK_ThrowRequireMsg(partition->num_entities() == 0, "Partition is not empty.");

    partition->m_partitionId = INVALID_PARTITION_ID;
    partition->reset_partition_id_in_owned_buckets();
    m_repoMeta.update_repo_meta_partition_removed(*partition);
    set_need_sync_from_partition(partition->get_rank(), true);
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
      set_need_sync_from_partition(rank, false);
    }

    sync_num_buckets_host_to_device();
    sync_num_partitions_host_to_device();
  }

  void sync_from_partitions(EntityRank rank)
  {
    Kokkos::Profiling::pushRegion("sync_from_partiitons");
    if (need_sync_from_partition(rank)) {
      sort_entities_in_buckets(rank);

      sort_buckets_in_partitions(rank);

      sort_partitions_and_sync_ids(rank);

      construct_bucket_proxy_list_and_copy_from_partitions(rank);

      // FIXME this step could be postponed and done later if preferred.
      copy_bucket_proxy_view_to_bucket_view(rank);

      update_fast_mesh_indices(rank);
    }
    Kokkos::Profiling::popRegion();

    // TODO
    // if (m_mesh->get_bulk_on_host().should_sort_buckets_by_first_entity_identifier()) {
    // }
  }

  // TODO: convert this to nested parallelism and have buckets' sortings
  //       be done in parallel
  void sort_entities_in_buckets(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto numBuckets = get_active_bucket_span(rank);

    for (unsigned i = 0; i < numBuckets; ++i) {
      auto& bucket = buckets(i);
      if (!bucket.is_active()) { continue; }

      if (bucket.size() == 0) {
        invalidate_bucket(&bucket);
      } else {
        bucket.sort_entities();
      }
    }
    Kokkos::fence();
  }

  // TODO: convert this to nested parallelism and have partitions' sortings
  //       be done in parallel
  void sort_buckets_in_partitions(EntityRank rank)
  {
    auto& partitions = m_partitions[rank];

    for (unsigned i = 0; i < get_active_partition_span(rank); ++i) {
      auto& partition = partitions(i);
      if (!partition.is_active()) { continue; }

      if (partition.num_buckets() == 0) {
        invalidate_partition(&partition);
      } else {
        partition.sort_buckets();
      }
    }
  }

  void sort_partitions_and_sync_ids(EntityRank rank)
  {
    auto& partitions = m_partitions[rank];
    DevicePartitionUView compactPartitionUView(partitions.data(), get_active_partition_span(rank));

    if (!Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactPartitionUView)) {
      Kokkos::sort(compactPartitionUView);
    }
    
    m_partitions[rank] = partitions;
    m_repoMeta.set_to_match_active_partitions(rank);

    for (unsigned i = 0; i < num_partitions(rank); ++i) {
      auto& partition = m_partitions[rank](i);
      partition.set_partition_id(i);
      partition.reset_partition_id_in_owned_buckets();
    }
  }

  void construct_bucket_proxy_list_and_copy_from_partitions(EntityRank rank)
  {
    auto numBuckets = num_buckets(rank);
    if (m_bucketsProxy[rank].extent(0) != numBuckets) {
      m_bucketsProxy[rank] = DeviceBucketWrapperView(Kokkos::view_alloc("DeviceBucketWrapperView", Kokkos::WithoutInitializing), numBuckets);
    }

    auto numPartitions = num_partitions(rank);
    auto& partitions = get_partitions(rank);
    auto& bucketsProxy = m_bucketsProxy[rank];

    for (unsigned i = 0, bucketIdx = 0; i < numPartitions; ++i) {
      auto& partition = partitions(i);
      auto numBucketInPartition = partition.num_buckets();
      for (unsigned j = 0; j < numBucketInPartition; ++j) {
        bucketsProxy(bucketIdx) = partition.m_buckets(j);
        bucketsProxy(bucketIdx)->m_bucketId = bucketIdx;
        bucketIdx++;
      }
    }
  }

  void update_fast_mesh_indices(EntityRank rank)
  {
    auto& fastMeshIndex = m_mesh->get_fast_mesh_indices();
    auto numBuckets = num_buckets(rank);
    auto& buckets = m_buckets[rank];

    Kokkos::parallel_for(numBuckets,
      KOKKOS_LAMBDA(const int bucketIdx) {
        auto& bucket = buckets(bucketIdx);

        if (!bucket.is_active()) { return; }

        for (unsigned i = 0; i < bucket.size(); ++i) {
          auto entity = bucket[i];

          if (!entity.is_local_offset_valid()) { continue; }
          fastMeshIndex(entity.local_offset()) = FastMeshIndex{bucket.bucket_id(), i};
        }
      }
    );

    Kokkos::fence();
  }

  void copy_bucket_proxy_view_to_bucket_view(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto& proxyBuckets = m_bucketsProxy[rank];
    auto numBuckets = num_buckets(rank);
    auto initCapacity = std::max(numBuckets, initialDeviceBucketViewCapacity);
    DeviceBucketView copyBucketView = DeviceBucketView(Kokkos::view_alloc("DeviceBucketView", Kokkos::WithoutInitializing), initCapacity);

    Kokkos::parallel_for(numBuckets,
      KOKKOS_LAMBDA(const int idx) {
        copyBucketView(idx) = *(proxyBuckets(idx).bucketPtr);
      }
    );
    Kokkos::fence();
    m_buckets[rank] = copyBucketView;

    m_repoMeta.set_to_match_active_buckets(rank);
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

  KOKKOS_INLINE_FUNCTION
  unsigned get_bucket_capacity() const { return m_maximumBucketCapacity; }
  
  KOKKOS_INLINE_FUNCTION
  unsigned get_initial_bucket_capacity() const { return m_initialBucketCapacity; }

  KOKKOS_INLINE_FUNCTION
  unsigned get_maximum_bucket_capacity() const { return m_maximumBucketCapacity; }

  KOKKOS_INLINE_FUNCTION
  unsigned num_buckets(EntityRank rank) const {
    return m_repoMeta.num_buckets(rank);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned num_partitions(EntityRank rank) const {
    return m_repoMeta.num_partitions(rank);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_last_active_bucket_idx(EntityRank rank) const {
    return m_repoMeta.get_last_active_bucket_idx(rank);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_last_active_partition_idx(EntityRank rank) const {
    return m_repoMeta.get_last_active_partition_idx(rank);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_active_bucket_span(EntityRank rank) const {
    return m_repoMeta.get_active_bucket_span(rank);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_active_partition_span(EntityRank rank) const {
    return m_repoMeta.get_active_partition_span(rank);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_next_avail_bucket_idx(EntityRank rank) const
  {
    return m_repoMeta.get_next_avail_bucket_idx(rank);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_next_avail_partition_idx(EntityRank rank) const {
    return m_repoMeta.get_next_avail_partition_idx(rank);
  }

  KOKKOS_INLINE_FUNCTION
  bool need_sync_from_partition(EntityRank rank) const {
    return m_repoMeta.m_needSyncFromPartitions[rank];
  }

  KOKKOS_INLINE_FUNCTION
  void set_need_sync_from_partition(EntityRank rank, bool needSync) {
    m_repoMeta.set_need_sync_from_partition(rank, needSync);
  }

  KOKKOS_INLINE_FUNCTION
  bool is_last_partition_reference(unsigned rank = stk::topology::NODE_RANK) const {
    return (m_partitions[rank].use_count() == 1);
  }

  KOKKOS_INLINE_FUNCTION
  bool is_last_bucket_reference(unsigned rank = stk::topology::NODE_RANK) const {
    return (m_buckets[rank].use_count() == 1);
  }

  KOKKOS_INLINE_FUNCTION
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

  template<typename EntityViewType>
  unsigned get_max_num_parts_per_entity(EntityViewType const& entities) const;

  unsigned m_initialBucketCapacity;
  unsigned m_maximumBucketCapacity;

  DeviceMeshT<NgpMemSpace>* m_mesh;
  DeviceBucketRepoMetadata<NgpMemSpace> m_repoMeta;
  
  DeviceBucketView m_buckets[stk::topology::NUM_RANKS];
  DevicePartitionView m_partitions[stk::topology::NUM_RANKS];

  // Used for intermediate bucket sorting without copying buckets
  DeviceBucketWrapperView m_bucketsProxy[stk::topology::NUM_RANKS];
};

template <typename NgpMemSpace>
template <typename EntityViewType, typename PartOrdinalProxyIndicesViewType>
void DeviceBucketRepository<NgpMemSpace>::batch_get_partitions(EntityViewType const& entities, PartOrdinalProxyIndicesViewType const& partOrdinalProxy, UnsignedPairViewType<NgpMemSpace>& srcDestPartitionIds) const
{
  auto& deviceMesh = *m_mesh;
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_CLASS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto rank = deviceMesh.entity_rank(entity);
      auto fastMeshIndex = deviceMesh.device_mesh_index(entity);
      auto bucketId = fastMeshIndex.bucket_id;
      auto& bucket = deviceMesh.get_bucket(rank, bucketId);
      auto srcPartitionId = bucket.m_owningPartitionId;
      auto& partitions = m_partitions[rank];

      auto span = get_active_partition_span(rank);
      srcDestPartitionIds(idx).first = srcPartitionId;

      // TODO convert to a proper search after tuning into a nested team-based loop
      for (unsigned i = 0; i < span; ++i) {
        auto& partition = partitions(i);
        auto partitionId = partition.partition_id();
        if (partitionId == INVALID_PARTITION_ID) { continue; }

        auto partitionPartOrdinals = partition.superset_part_ordinals();
        if (partOrdinalProxy(idx).length != partitionPartOrdinals.extent(0)) { continue; }

        auto foundId = partitionId;
        for (unsigned j = 0; j < partOrdinalProxy(idx).length; j++) {
          auto localPartOrdinal = partOrdinalProxy(idx).startPtr + j;

          if (*localPartOrdinal != partitionPartOrdinals(j)) {
            foundId = INVALID_PARTITION_ID;
            break;
          }
        }
        if (foundId != INVALID_PARTITION_ID) {
          srcDestPartitionIds(idx).second = foundId;
          break;
        }
      }
    }
  );

  Kokkos::fence();
}

} } }

#endif
