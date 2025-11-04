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
#include "DevicePartition.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Macros.hpp"
#include "NgpMeshImpl.hpp"
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
#include "stk_mesh/baseImpl/DeviceMeshViewVector.hpp"
#include "stk_mesh/baseImpl/DevicePartition.hpp"
#include "stk_mesh/baseImpl/NgpMeshImpl.hpp"
#include "stk_mesh/baseImpl/ViewVector.hpp"
#include "Kokkos_Core.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace mesh {
  
template <typename NgpMemSpace>
class DeviceMeshT;

namespace impl {

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
  using DeviceBucketViewVector = ImplDeviceMeshViewVector<DeviceBucket, stk::ngp::UVMMemSpace>;
  using DeviceBucketUView = Kokkos::View<DeviceBucket*, stk::ngp::UVMMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using DevicePartitionViewVector = ImplDeviceMeshViewVector<DevicePartition<NgpMemSpace>, stk::ngp::UVMMemSpace>;
  using DevicePartitionUView = Kokkos::View<DevicePartition<NgpMemSpace>*, stk::ngp::UVMMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using DevicePartRankView = EntityRankViewType<NgpMemSpace>;
  using BoolViewType  = Kokkos::View<bool*, stk::ngp::MemSpace>;

  KOKKOS_DEFAULTED_FUNCTION
  DeviceBucketRepository() = default;

  DeviceBucketRepository(DeviceMeshT<NgpMemSpace>* deviceMesh,
                         unsigned initialBucketCapacity = get_default_initial_bucket_capacity(),
                         unsigned maximumBucketCapacity = get_default_maximum_bucket_capacity())
    : m_initialBucketCapacity(initialBucketCapacity),
      m_maximumBucketCapacity(maximumBucketCapacity),
      m_mesh(deviceMesh),
      m_needSyncFromPartitions{},
      m_numInternalParts(0),
      m_lastInternalPartOrdinal(InvalidPartOrdinal)
  {}

  KOKKOS_DEFAULTED_FUNCTION DeviceBucketRepository(DeviceBucketRepository const&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceBucketRepository(DeviceBucketRepository&&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceBucketRepository& operator=(DeviceBucketRepository&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceBucketRepository& operator=(DeviceBucketRepository const&) = default;
  KOKKOS_DEFAULTED_FUNCTION ~DeviceBucketRepository() = default;

  KOKKOS_INLINE_FUNCTION
  DevicePartitionViewVector& get_partitions(const EntityRank rank)
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
      STK_ThrowRequireMsg(partitionId < m_partitions[rank].size(), "Invalid partition id.");
    ))

    return &(m_partitions[rank][partitionId]);
  }
      
  DevicePartition<NgpMemSpace>* get_partition(const EntityRank rank,
                                              PartOrdinalViewType<NgpMemSpace> const& partOrdinals) const
  {
    Kokkos::Profiling::pushRegion("get_partition");
    STK_ThrowRequireMsg(rank <= stk::topology::NUM_RANKS, "Invalid rank = " << rank);

    auto policy = Kokkos::RangePolicy(0, m_partitions[rank].size());
    DevicePartitionUView compactPartitionUView(m_partitions[rank].data(), m_partitions[rank].size());
    auto partitionId = search_matching_device_partitions(policy, compactPartitionUView, partOrdinals);
    Kokkos::Profiling::popRegion();

    return (partitionId != INVALID_BUCKET_ID) ? &(m_partitions[rank][partitionId]) : nullptr;
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

    auto& partitions = m_partitions[rank];
    auto newPartitionIdx = partitions.size();
    partitions.emplace_back(m_mesh, this, partOrdinals, rank, newPartitionIdx, copyPartOrdinals);
    Kokkos::Profiling::popRegion();

    return &partitions[partitions.size()-1];
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
        if (srcPartitionId == srcDestPartitions(i).second) {
          return;
        }

        auto& partition = m_partitions[rank][srcPartitionId];
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
      if (destPartitionId == hostSrcDestPartitions(i).first) {
        continue;
      }

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
      STK_ThrowRequireMsg(static_cast<unsigned>(bucketIdx) < m_buckets[rank].size(), "Invalid bucket id.");
    ))

    return &(m_buckets[rank][bucketIdx]);
  }

  DeviceBucket* construct_new_bucket_without_part_ordinals(const EntityRank rank)
  {
    DeviceBucketViewVector& buckets = m_buckets[rank];

    // if bucket view vector is resized, bucket ptrs in partitions need to be reassigned
    bool needToResize = buckets.size() == buckets.capacity();

    buckets.emplace_back();
    DeviceBucket& newBucket = buckets[buckets.size()-1];

    newBucket.m_bucketSize = 0;
    newBucket.m_bucketCapacity = get_bucket_capacity();
    newBucket.m_bucketId = buckets.size()-1;
    newBucket.m_entityRank = rank;

    if (needToResize) {
      update_bucket_ptrs_in_partitions(rank);
    }
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
    Kokkos::fence();
  }

  void invalidate_bucket(DeviceBucket* bucket)
  {
    STK_ThrowRequireMsg(bucket != nullptr, "Invalid bucket pointer.");
    STK_ThrowRequireMsg(bucket->size() == 0, "Bucket is not empty.");
  
    bucket->m_bucketId = INVALID_BUCKET_ID;
    auto partition = get_partition(bucket->entity_rank(), bucket->partition_id());
    partition->update_partition_meta_bucket_removed(bucket);
    m_buckets[bucket->entity_rank()].decrement_num_active_entries();
    set_need_sync_from_partition(bucket->entity_rank(), true);
  }

  void invalidate_partition(DevicePartition<NgpMemSpace>* partition)
  {
    STK_ThrowRequireMsg(partition != nullptr, "Invalid partition pointer.");
    STK_ThrowRequireMsg(partition->num_buckets() == 0, "Partition is not empty.");
    STK_ThrowRequireMsg(partition->num_entities() == 0, "Partition is not empty.");

    partition->m_partitionId = INVALID_PARTITION_ID;
    partition->reset_partition_id_in_owned_buckets();
    m_partitions[partition->get_rank()].decrement_num_active_entries();
    set_need_sync_from_partition(partition->get_rank(), true);
  }

  void sync_from_partitions()
  {
    for (auto rank = stk::topology::NODE_RANK; rank < stk::topology::NUM_RANKS; ++rank) {
      sync_from_partitions(rank);
      clear_need_sync_from_partition(rank);
    }
  }

  void sync_from_partitions(EntityRank rank)
  {
    Kokkos::Profiling::pushRegion("sync_from_partiitons");
    if (need_sync_from_partition(rank)) {
      sort_entities_in_buckets(rank);

      sort_buckets_in_partitions(rank);

      sort_partitions_and_sync_ids(rank);

      copy_sorted_buckets_from_partitions(rank);

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
    auto numBuckets = buckets.size();

    for (unsigned i = 0; i < numBuckets; ++i) {
      auto& bucket = buckets[i];
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
    auto span = partitions.size();

    for (unsigned i = 0; i < span; ++i) {
      auto& partition = partitions[i];
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
    Kokkos::Profiling::pushRegion("sort_partitions_and_sync_ids");
    auto& partitions = get_partitions(rank);
    DevicePartitionUView compactPartitionUView(partitions.data(), partitions.size());

    bool isSorted = Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactPartitionUView);
    if (!isSorted) {
      Kokkos::sort(compactPartitionUView);
    }

    STK_ThrowAssert(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactPartitionUView));

    for (unsigned i = 0; i < num_partitions(rank); ++i) {
      auto& partition = m_partitions[rank][i];
      partition.set_partition_id(i);
      partition.reset_partition_id_in_owned_buckets();
    }
  }

  void copy_sorted_buckets_from_partitions(EntityRank rank)
  {
    Kokkos::Profiling::pushRegion("sort_buckets");
    auto bucketSpan = m_buckets[rank].size();
    DeviceBucketViewVector copyBucketView("DeviceBucketView", bucketSpan);
    DeviceBucketUView uview(m_buckets[rank].data(), bucketSpan);

    auto execSpace = typename NgpMemSpace::execution_space{};
    Kokkos::Experimental::move(execSpace, uview, copyBucketView.get_view());

    auto numPartitions = num_partitions(rank);
    auto& partitions = get_partitions(rank);
    auto& buckets = m_buckets[rank];

    for (unsigned i = 0, bucketIdx = 0; i < numPartitions; ++i) {
      auto& partition = partitions[i];

      if (!partition.is_active()) { continue; }

      auto numBucketsInPartition = partition.num_buckets();
      for (unsigned j = 0; j < numBucketsInPartition; ++j) {
        auto bucketIdFromPartition = partition.m_buckets[j].bucketId;
        auto& bucketToCopy = copyBucketView[bucketIdFromPartition];

        STK_ThrowAssertMsg(bucketIdFromPartition == bucketToCopy.bucket_id(), "Mismatched copying bucket index and bucket id");
        buckets[bucketIdx] = bucketToCopy;
        buckets[bucketIdx].set_bucket_id(bucketIdx);
        partition.m_buckets[j].bucketPtr = &buckets[bucketIdx];
        partition.m_buckets[j].bucketId = bucketIdx;
        bucketIdx++;
      }
    }

    Kokkos::Profiling::popRegion();
  }

  void update_fast_mesh_indices(EntityRank rank)
  {
    auto& fastMeshIndex = m_mesh->get_fast_mesh_indices();
    auto numBuckets = num_buckets(rank);
    auto& buckets = m_buckets[rank];

    Kokkos::parallel_for(numBuckets,
      KOKKOS_LAMBDA(const int bucketIdx) {
        auto& bucket = buckets[bucketIdx];
        auto bucketSize = bucket.size();

        if (!bucket.is_active()) { return; }

        for (unsigned i = 0; i < bucketSize; ++i) {
          auto entity = bucket[i];
          if (!entity.is_local_offset_valid()) { continue; }
          fastMeshIndex(entity.local_offset()) = FastMeshIndex{bucket.bucket_id(), i};
        }
      }
    );

    Kokkos::fence();
  }

  void reset_and_invalidate_all_buckets(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto bucketCount = buckets.size();

    for (unsigned i = 0; i < bucketCount; ++i) {
      auto& bucket = buckets[i];
      bucket.m_bucketId = INVALID_BUCKET_ID;
    }
    buckets.clear_active_entries();
  }

  void reset_and_invalidate_all_partitions(EntityRank rank)
  {
    auto& partitions = m_partitions[rank];
    auto partitionCount = partitions.size();

    for (unsigned i = 0; i < partitionCount; ++i) {
      auto& partition = partitions[i];
      partition.m_partitionId = INVALID_PARTITION_ID;
    }
    partitions.clear_active_entries();
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
    return m_buckets[rank].num_active_entries();
  }

  KOKKOS_INLINE_FUNCTION
  unsigned num_partitions(EntityRank rank) const {
    return m_partitions[rank].num_active_entries();
  }

  bool need_sync_from_partition(EntityRank rank) const {
    return m_needSyncFromPartitions[rank];
  }

  void set_need_sync_from_partition(EntityRank rank, bool needSync) {
    m_needSyncFromPartitions[rank] = m_needSyncFromPartitions[rank] | needSync;
  }

  void clear_need_sync_from_partition(EntityRank rank) {
    m_needSyncFromPartitions[rank] = false;
  }

  KOKKOS_INLINE_FUNCTION
  bool is_last_partition_reference(unsigned rank = stk::topology::NODE_RANK) const {
    return (m_partitions[rank].get_view().use_count() == 1);
  }

  KOKKOS_INLINE_FUNCTION
  bool is_last_bucket_reference(unsigned rank = stk::topology::NODE_RANK) const {
    return (m_buckets[rank].get_view().use_count() == 1);
  }

  KOKKOS_INLINE_FUNCTION
  void clear_device_buckets_and_partitions()
  {
    KOKKOS_IF_ON_HOST((
      for (EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::END_RANK; rank++) {
        if (is_last_partition_reference(rank)) {
          auto partitions = std::max(m_partitions[rank].size(), num_partitions(rank));
          for (unsigned iPartition = 0; iPartition < partitions; ++iPartition) {
            m_partitions[rank][iPartition].~DevicePartition();
          }
          if (is_last_bucket_reference(rank)) {
            auto buckets = std::max(m_buckets[rank].size(), num_buckets(rank));
            for (unsigned iBucket = 0; iBucket < buckets; ++iBucket) {
              m_buckets[rank][iBucket].~DeviceBucket();
            }
          }
        }
      }
    ))
  }

  KOKKOS_INLINE_FUNCTION
  EntityRank get_part_rank(unsigned partOrdinal) const
  {
    KOKKOS_IF_ON_HOST((
      auto& hostPart = m_mesh->get_bulk_on_host().mesh_meta_data().get_part(partOrdinal);
      return hostPart.primary_entity_rank();
    ))
    KOKKOS_IF_ON_DEVICE((
      return m_partRanks(partOrdinal);
    ))
  }

  KOKKOS_INLINE_FUNCTION
  bool is_ranked_part(unsigned partOrdinal) const
  {
    return get_part_rank(partOrdinal) != stk::topology::INVALID_RANK;
  }

  KOKKOS_INLINE_FUNCTION
  unsigned is_internal_part(unsigned partOrdinal) const
  {
    if (partOrdinal == InvalidPartOrdinal || m_lastInternalPartOrdinal == InvalidPartOrdinal) {
      Kokkos::abort("Part ordinal must not be invalid.\n");
    }

    return partOrdinal <= m_lastInternalPartOrdinal;
  }

  void copy_device_part_rank_info_from_host()
  {
    auto& allParts = m_mesh->get_bulk_on_host().mesh_meta_data().get_parts();
    Kokkos::realloc(m_partRanks, allParts.size());

    auto hostPartRanks = Kokkos::create_mirror_view(m_partRanks);

    for (unsigned i = 0; i < allParts.size(); ++i) {
      auto& part = *allParts[i];

      if (stk::mesh::impl::is_internal_part(part)) {
        m_numInternalParts++;
      }

      hostPartRanks(i) = part.primary_entity_rank();
    }
    Kokkos::deep_copy(m_partRanks, hostPartRanks);
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
    DeviceBucketViewVector bucketBuffer("BucketBuffer");
    bucketBuffer.reserve(bucketCapacityToInit);

    auto partitionCapacityToInit = hostPartitions.size();
    DevicePartitionViewVector partitionBuffer("PartitionBuffer");
    partitionBuffer.reserve(partitionCapacityToInit);

    for (unsigned partitionIdx = 0, bucketIdxInView = 0; partitionIdx < numHostPartitions; ++partitionIdx) {
      auto hostPartition = hostPartitions[partitionIdx];
      auto& partitionParts = hostPartition->get_legacy_partition_id();
      auto numBucketsInPartition = hostPartition->num_buckets();
      partitionBuffer.emplace_back(m_mesh, this, partitionParts, rank, partitionIdx);

      for (unsigned i = 0; i < numBucketsInPartition; ++i, ++bucketIdxInView) {
        auto& iBucket = bucketIdxInView;
        stk::mesh::Bucket& stkBucket = *hostBuckets[iBucket];
        const unsigned ngpBucketId = stkBucket.ngp_mesh_bucket_id();

        if (ngpBucketId == INVALID_BUCKET_ID) {
          bucketBuffer.emplace_back();

          bucketBuffer[iBucket].initialize_bucket_attributes(stkBucket);
          bucketBuffer[iBucket].initialize_part_ordinals_from_host(stkBucket);
          bucketBuffer[iBucket].update_entity_data_from_host(stkBucket);
          bucketBuffer[iBucket].update_sparse_connectivity_from_host(stkBucket);
          anyBucketChanges = true;
        } else {
          bucketBuffer.emplace_back(m_buckets[rank][ngpBucketId]);

          if (stkBucket.ngp_mesh_bucket_is_modified()) {
            bucketBuffer[iBucket].update_entity_data_from_host(stkBucket);
            bucketBuffer[iBucket].update_sparse_connectivity_from_host(stkBucket);
            anyBucketChanges = true;
          }
          bucketBuffer[iBucket].m_bucketId = stkBucket.bucket_id();
        }

        stkBucket.set_ngp_mesh_bucket_id(iBucket);

        // Update device buckets in partition
        partitionBuffer[partitionIdx].add_bucket(&bucketBuffer[iBucket]);

        STK_ThrowRequire(i == partitionBuffer[partitionIdx].num_buckets()-1);
        STK_ThrowRequire(iBucket == stkBucket.bucket_id());
      }
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

    copy_device_part_rank_info_from_host();
    Kokkos::Profiling::popRegion();
  }

  template<typename EntityViewType>
  unsigned get_max_num_parts_per_entity(EntityViewType const& entities) const;

  void update_bucket_ptrs_in_partitions(EntityRank rank)
  {
    auto& buckets = m_buckets[rank];
    auto& partitions = m_partitions[rank];

    auto numActiveBucketSpan = buckets.size();
    Kokkos::parallel_for(numActiveBucketSpan,
      KOKKOS_CLASS_LAMBDA(const int idx) {
        auto& bucket = buckets[idx];

        if (!bucket.is_active()) { return; }

        auto partitionId = bucket.partition_id();

        // Newly created bucket
        if (partitionId == INVALID_PARTITION_ID) { return; }

        auto& partition = partitions[partitionId];
        partition.update_bucket_ptrs(bucket);
      }
    );
    Kokkos::fence();
  }

  void update_last_internal_part_ordinal()
  {
    auto& meta = m_mesh->get_bulk_on_host().mesh_meta_data();
    auto& allParts = meta.get_parts();
    bool updated = false;

    for (auto part : allParts) {
      if (part != nullptr) {
        if (stk::mesh::impl::is_internal_part(*part)) {
          m_lastInternalPartOrdinal = part->mesh_meta_data_ordinal();
          updated = true;
        }
      }
    }

    if (!updated)
      m_lastInternalPartOrdinal = InvalidPartOrdinal;
  }

  unsigned m_initialBucketCapacity;
  unsigned m_maximumBucketCapacity;

  DeviceMeshT<NgpMemSpace>* m_mesh;

  DeviceBucketViewVector m_buckets[stk::topology::NUM_RANKS];
  DevicePartitionViewVector m_partitions[stk::topology::NUM_RANKS];
  bool m_needSyncFromPartitions[stk::topology::NUM_RANKS];

  // TODO Refactor to use a proper device side part (either DevicePart or a device copyable Part)
  DevicePartRankView m_partRanks;
  unsigned m_numInternalParts;
  unsigned m_lastInternalPartOrdinal;
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

      auto span = m_partitions[rank].size();
      srcDestPartitionIds(idx).first = srcPartitionId;

      // TODO convert to a proper search after tuning into a nested team-based loop
      for (unsigned i = 0; i < span; ++i) {
        auto& partition = partitions[i];
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
