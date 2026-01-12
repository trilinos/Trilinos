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

struct EntitySrcDest {
  Entity entity;
  EntityRank rank;

  unsigned srcPartitionId;
  unsigned srcBucketId;
  unsigned srcBucketOrd;

  unsigned destPartitionId;
  unsigned destBucketIndexInPartition;
  unsigned destBucketOrd;

  unsigned destBucketId = INVALID_BUCKET_ID;
  unsigned destNodeConnectivityStartIdx = INVALID_INDEX;
  unsigned destNodeConnectivityLength = 0 ;
  unsigned destSparseConnectivityStartIdx = INVALID_INDEX;
  unsigned destSparseConnectivityLength = 0;

  KOKKOS_INLINE_FUNCTION
  bool operator<(EntitySrcDest const& rhs) const {
    if (rank != rhs.rank) {
      return rank < rhs.rank;
    } else {
     return destPartitionId < rhs.destPartitionId;
    }
  }
};

struct BucketConnectivitySizes {
  unsigned bucketId;
  unsigned partitionId;
  unsigned bucketIndexInPartition;

  unsigned totalNumEntities;
  unsigned newNodeConnectivityViewSize;
  unsigned newSparseConnectivityViewSize;
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
      endRank(deviceMesh->get_end_rank()),
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

  template <typename EntityViewType, typename PartOrdinalProxyIndicesViewType, typename EntitySrcDestView>
  void batch_get_partitions(EntityViewType const& entities,
      PartOrdinalProxyIndicesViewType const& partOrdinalProxy,
      EntitySrcDestView& entitySrcDestView) const;

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

  template <typename EntitySrcDestView>
  void batch_move_entities(EntitySrcDestView const& entitySrcDestView)
  {
    // TODO: convert this to a nested parallel for with a team policy
    // so that copy_bucket_connectivity_from_src can use parallelized copy
    Kokkos::parallel_for(
        entitySrcDestView.extent(0), KOKKOS_CLASS_LAMBDA(const int i) {
          auto& entitySrcDest = entitySrcDestView(i);
          auto& entity = entitySrcDest.entity;
          auto srcPartitionId = entitySrcDest.srcPartitionId;
          auto srcBucketId = entitySrcDest.srcBucketId;
          auto srcBucketOrd = entitySrcDest.srcBucketOrd;
          auto destPartitionId = entitySrcDest.destPartitionId;
          auto destBucketIndexInPartition = entitySrcDest.destBucketIndexInPartition;
          auto destBucketId = entitySrcDest.destBucketId;
          auto destBucketOrd = entitySrcDest.destBucketOrd;
          auto rank = entitySrcDest.rank;

          auto srcPartition = get_partition(rank, srcPartitionId);
          auto destPartition = get_partition(rank, destPartitionId);
          auto srcBucket = get_bucket(rank, srcBucketId);
          auto destBucket = get_bucket(rank, destBucketId);

          destPartition->add_entity_without_connectivity(entity, destBucketIndexInPartition, destBucketOrd);

          copy_bucket_connectivities_from_src(srcBucket, destBucket, entitySrcDest, endRank);
          srcPartition->remove_entity(entity, srcBucketId, srcBucketOrd);

          srcPartition->set_modified();
          destPartition->set_modified();
        });

    Kokkos::fence();

    for(EntityRank entityRank = stk::topology::BEGIN_RANK; entityRank < endRank; ++entityRank) {
      move_field_data(entityRank);
    }
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

  template <typename NewBucketsToAddViewType>
  void batch_create_buckets(NewBucketsToAddViewType numNewBucketsToAddInPartition)
  {
    auto hostNumBucketsToCreate =
        Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, numNewBucketsToAddInPartition);
    for (unsigned i = 0; i < hostNumBucketsToCreate.extent(0); ++i) {
      auto& numBucktsToCreateInfo = hostNumBucketsToCreate(i);
      auto rank = numBucktsToCreateInfo.rank;
      auto partitionId = numBucktsToCreateInfo.partitionId;
      auto numBucketsToCreate = numBucktsToCreateInfo.numBucketsToAdd;
      if (numBucketsToCreate == 0 || partitionId == INVALID_PARTITION_ID || rank == topology::INVALID_RANK) {
        continue;
      }

      auto partition = get_partition(rank, partitionId);
      for (unsigned j = 0; j < numBucketsToCreate; ++j) {
        partition->create_new_bucket_without_connectivities();
        partition->set_modified();
      }
    }
  }

  template <typename EntitySrcDestView, typename BucketConnectivitySize2DView>
  void batch_init_bucket_connectivity_views(EntitySrcDestView const& entitySrcDestView,
                                            BucketConnectivitySize2DView& numNewConnectivityViewSizesInBuckets)
  {
    auto hostEntitySrcDestView = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, entitySrcDestView);
    auto hostNumNewConnectivitySizes = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, numNewConnectivityViewSizesInBuckets);

    for (unsigned i = 0; i < hostEntitySrcDestView.extent(0); ++i) {
      auto& entitySrcDest = hostEntitySrcDestView(i);
      auto entity = entitySrcDest.entity;
      auto rank = entitySrcDest.rank;
      auto bucketId = entitySrcDest.destBucketId;
      auto bucket = get_bucket(rank, bucketId);
      auto& connectivitySizes = hostNumNewConnectivitySizes(rank, bucketId);

      Kokkos::resize(Kokkos::WithoutInitializing, bucket->m_nodeConnectivityOffsets, connectivitySizes.totalNumEntities+1);
      Kokkos::resize(Kokkos::WithoutInitializing, bucket->m_nodeConnectivity, connectivitySizes.newNodeConnectivityViewSize);

      auto numNodes = m_mesh->get_bulk_on_host().num_nodes(entity);
      if (numNodes > bucket->m_nodeOrdinals.extent(0)) {
        Kokkos::resize(Kokkos::WithoutInitializing, bucket->m_nodeOrdinals, numNodes);

        Kokkos::parallel_for(numNodes, KOKKOS_LAMBDA(const int idx) { bucket->m_nodeOrdinals(idx) = idx; });
      }

      Kokkos::resize(Kokkos::WithoutInitializing, bucket->m_sparseConnectivityOffsets, stk::topology::NUM_RANKS+1, connectivitySizes.totalNumEntities);
      Kokkos::resize(Kokkos::WithoutInitializing, bucket->m_sparseConnectivity, connectivitySizes.newSparseConnectivityViewSize);
      Kokkos::resize(Kokkos::WithoutInitializing, bucket->m_sparseConnectivityOrdinals, connectivitySizes.newSparseConnectivityViewSize);

      auto srcBucket = get_bucket(rank, entitySrcDest.srcBucketId);
      if (srcBucket->has_permutations()) {
        bucket->set_has_permutations(true);
        Kokkos::resize(Kokkos::WithoutInitializing, bucket->m_sparseConnectivityPermutations, connectivitySizes.newSparseConnectivityViewSize);
      }
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

  // TODO refactor
  DeviceBucket* construct_new_bucket_without_part_ordinals(const EntityRank rank)
  {
    DeviceBucketViewVector& buckets = m_buckets[rank];

    // if bucket view vector is resized, bucket ptrs in partitions need to be reassigned
    bool needToResize = buckets.size() == buckets.capacity();

    buckets.emplace_back();
    DeviceBucket& newBucket = buckets[buckets.size()-1];

    newBucket.m_owningMesh = m_mesh;
    newBucket.m_bucketId = buckets.size()-1;
    newBucket.m_bucketSize = 0;
    newBucket.m_bucketCapacity = get_bucket_capacity();
    newBucket.m_entityRank = rank;
    newBucket.m_endRank = endRank;

    if (needToResize) {
      update_bucket_ptrs_in_partitions(rank);
    }
    newBucket.set_modified();

    return &newBucket;
  }

  DeviceBucket* construct_new_bucket(const EntityRank rank,
                                     PartOrdinalViewType<NgpMemSpace> const& partOrdinals)
  {
    Kokkos::Profiling::pushRegion("construct_bucket");
    auto newBucket = construct_new_bucket_without_part_ordinals(rank);
    newBucket->m_partOrdinals = partOrdinals;
    DeviceFieldDataManagerBase* deviceFieldDataManagerBase = m_mesh->get_field_data_manager(m_mesh->get_bulk_on_host());

    auto hostPartOrdinals = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, partOrdinals);

    const MetaData& meta = m_mesh->get_bulk_on_host().mesh_meta_data();
    const PartVector& allParts = meta.get_parts();
    PartVector parts(partOrdinals.extent(0));
    for(unsigned i=0; i<partOrdinals.extent(0); ++i) {
      parts[i] = allParts[hostPartOrdinals(i)];
    }

    constexpr bool deviceMeshMod = true;
    deviceFieldDataManagerBase->add_new_bucket(rank, newBucket->size(), newBucket->capacity(), parts, deviceMeshMod);

    auto& fieldsOfRank = meta.get_fields(rank);
    for(auto* field : fieldsOfRank) {
      if (field->has_device_data()) {
        deviceFieldDataManagerBase->set_device_field_meta_data(*impl::get_device_data(*field));
      }
    }

    Kokkos::Profiling::popRegion();
    return newBucket;
  }

  void invalidate_bucket(DeviceBucket* bucket)
  {
    STK_ThrowRequireMsg(bucket != nullptr, "Invalid bucket pointer.");
    STK_ThrowRequireMsg(bucket->size() == 0, "Bucket is not empty.");

    bucket->m_bucketId = INVALID_BUCKET_ID;
    auto partition = get_partition(bucket->entity_rank(), bucket->partition_id());
    partition->update_partition_meta_bucket_removed(bucket);
    m_buckets[bucket->entity_rank()].decrement_num_active_entries();

    bucket->set_modified();
    partition->set_modified();
  }

  void invalidate_partition(DevicePartition<NgpMemSpace>* partition)
  {
    STK_ThrowRequireMsg(partition != nullptr, "Invalid partition pointer.");
    STK_ThrowRequireMsg(partition->num_buckets() == 0, "Partition is not empty.");
    STK_ThrowRequireMsg(partition->num_entities() == 0, "Partition is not empty.");

    partition->m_partitionId = INVALID_PARTITION_ID;
    partition->reset_partition_id_in_owned_buckets();
    m_partitions[partition->get_rank()].decrement_num_active_entries();

    partition->set_modified();
  }

  void sync_from_partitions()
  {
    Kokkos::Profiling::pushRegion("sync_from_partitions");
    auto needSyncFromPartitionPerRank = need_sync_from_partition();
    auto needSync = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, needSyncFromPartitionPerRank);

    for (auto rank = stk::topology::NODE_RANK; rank < stk::topology::NUM_RANKS; ++rank) {
      if (needSync(rank)) {
        sync_from_partitions(rank);
        clear_modified_from_all(rank);
      }
    }
    check_mesh_consistency();
    // print_part_ordinals();
    Kokkos::Profiling::popRegion();
  }

  void sync_from_partitions(EntityRank rank)
  {
    compress_connectivities_in_buckets(rank);

    sort_entities_in_buckets(rank);

    sort_buckets_in_partitions(rank);

    sort_partitions_and_sync_ids(rank);

    std::vector<BucketShift> bucketShifts = copy_sorted_buckets_from_partitions(rank);

    update_field_buckets(rank, bucketShifts);

    update_fast_mesh_indices(rank);

    // TODO
    // if (m_mesh->get_bulk_on_host().should_sort_buckets_by_first_entity_identifier()) {
    // }

  }


  void compress_connectivities_in_buckets(EntityRank rank);

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
    partitions.resize(partitions.num_active_entries());
    Kokkos::Profiling::popRegion();
  }

  using BucketShift = DeviceFieldDataManagerBase::BucketShift;

  std::vector<BucketShift> copy_sorted_buckets_from_partitions(EntityRank rank)
  {
    Kokkos::Profiling::pushRegion("sort_buckets");
    auto bucketSpan = m_buckets[rank].num_active_entries();
    DeviceBucketViewVector newBucketView("DeviceBucketView", bucketSpan);

    auto numPartitions = num_partitions(rank);
    auto& partitions = get_partitions(rank);
    auto& buckets = m_buckets[rank];
    std::vector<BucketShift> bucketShifts;
    for (unsigned i = 0, bucketIdx = 0; i < numPartitions; ++i) {
      auto& partition = partitions[i];

      if (!partition.is_active()) { continue; }

      auto numBucketsInPartition = partition.num_buckets();
      for (unsigned j = 0; j < numBucketsInPartition; ++j) {
        auto bucketIdFromPartition = partition.m_buckets[j].bucketId;

        STK_ThrowAssertMsg(bucketIdFromPartition == buckets[bucketIdFromPartition].bucket_id(), "Mismatched moving bucket index and bucket id");
        bucketShifts.emplace_back(bucketIdFromPartition, bucketIdx);
        newBucketView[bucketIdx] = buckets[bucketIdFromPartition];
        newBucketView[bucketIdx].set_bucket_id(bucketIdx);
        partition.m_buckets[j].bucketPtr = &newBucketView[bucketIdx];
        partition.m_buckets[j].bucketId = bucketIdx;
        bucketIdx++;
      }
    }
    newBucketView.set_active_entries(bucketSpan);
    m_buckets[rank] = newBucketView;

    buckets.resize(buckets.num_active_entries());

    Kokkos::Profiling::popRegion();

    return bucketShifts;
  }

  void move_field_data(EntityRank rank)
  {
    // Note: this only works for *moves* of entities, ie. each entity is removed from its
    //       old bucket and placed in a new one.  It does not work for more general
    //       sorting, ex. entity at index 1 is moved to index 2 and entity that was at
    //       index 2 is moved to index 3.  If the field data from entity at index 1 is copied
    //       to index 2, then the field data for the entity at index 2 is lost
    using FieldDataBytesType = FieldDataBytes<stk::ngp::DeviceSpace>;
    using FieldDataBytesViewType = Kokkos::View<FieldDataBytesType*, NgpMemSpace>;
    using InitValsPtrViewType = Kokkos::View<BytePtr*, stk::ngp::HostPinnedSpace>;

    FieldDataBytesViewType fieldDataBytesView(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fieldDataBytes"),0);

    InitValsPtrViewType initValsPtrsView;
    const MetaData& meta = m_mesh->get_bulk_on_host().mesh_meta_data();

    const FieldVector& rankFields = meta.get_fields(rank);
    if (!rankFields.empty()) {
      assemble_field_data_bytes_on_device<stk::ngp::DeviceSpace>(meta.get_fields(rank), fieldDataBytesView);
    }

    assemble_field_init_vals_on_device(rankFields, initValsPtrsView);

    const auto& fastMeshIndex = m_mesh->get_fast_mesh_indices();
    auto numBuckets = num_buckets(rank);
//    auto ngpMesh = *m_mesh;
    Kokkos::parallel_for(numBuckets,
      KOKKOS_CLASS_LAMBDA(const int bucketIdx) {
        auto bucket = get_bucket(rank, bucketIdx);
        auto bucketSize = bucket->size();

        if (!bucket->is_active()) { return; }

        for (unsigned i = 0; i < bucketSize; ++i) {
          auto entity = (*bucket)[i];
          if (!entity.is_local_offset_valid()) { continue; }

          FastMeshIndex oldFmi = fastMeshIndex(entity.local_offset());
          FastMeshIndex newFmi{bucket->bucket_id(), i};
          if (newFmi != oldFmi && fieldDataBytesView.extent(0) > 0) {
            copy_entity_bytes_kernel<stk::mesh::Layout::Left>(fieldDataBytesView, oldFmi, newFmi, initValsPtrsView);
          }

        }
      }
    );

    Kokkos::fence();
  }

  void update_fast_mesh_indices(EntityRank rank)
  {

    auto& fastMeshIndex = m_mesh->get_fast_mesh_indices();
    auto numBuckets = num_buckets(rank);
    Kokkos::parallel_for(numBuckets,
      KOKKOS_CLASS_LAMBDA(const int bucketIdx) {
        auto bucket = get_bucket(rank, bucketIdx);
        auto bucketSize = bucket->size();

        if (!bucket->is_active()) { return; }

        for (unsigned i = 0; i < bucketSize; ++i) {
          auto entity = (*bucket)[i];
          if (!entity.is_local_offset_valid()) { continue; }

          FastMeshIndex newFmi{bucket->bucket_id(), i};
          fastMeshIndex(entity.local_offset()) = newFmi;
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

  BoolViewType need_sync_from_partition() const
  {
    using ExecSpace = typename NgpMemSpace::execution_space;
    using TeamHandle = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

    BoolViewType needSyncPerRank("", stk::topology::NUM_RANKS);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<>(needSyncPerRank.extent(0), Kokkos::AUTO), KOKKOS_CLASS_LAMBDA(TeamHandle const& team) {
          auto rank = static_cast<stk::topology::rank_t>(team.league_rank());
          bool needSyncForThisRank = false;

          auto numPartitions = num_partitions(rank);
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, numPartitions),
              [&, this](const int i, bool& update) {
                auto partition = get_partition(rank, i);
                bool needSyncForThisPartition = partition->is_modified();

                if (!needSyncForThisPartition) {
                  for (unsigned j = 0; j < partition->num_buckets(); ++j) {
                    auto bucket = partition->m_buckets[j];
                    update |= bucket->is_modified();
                  }
                }
                update |= needSyncForThisPartition;
              },
              needSyncForThisRank);

          Kokkos::single(Kokkos::PerTeam(team), [&]() { needSyncPerRank(rank) = needSyncForThisRank; });
        });
    Kokkos::fence();

    return needSyncPerRank;
  }

  void clear_modified_from_all(const EntityRank rank)
  {
    using TeamHandle = Kokkos::TeamPolicy<>::member_type;

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<>(num_partitions(rank), Kokkos::AUTO), KOKKOS_CLASS_LAMBDA(TeamHandle const& team) {
          const int i = team.league_rank();
          auto partition = get_partition(rank, i);
          partition->clear_modified();

          auto numBuckets = partition->num_buckets();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, numBuckets), [&](const int j) {
            auto bucket = partition->m_buckets[j];
            bucket->clear_modified();
          });
        });
    Kokkos::fence();
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

  // FIXME basic way to check if an entity is valid
  KOKKOS_INLINE_FUNCTION
  bool is_valid(Entity const& entity) const
  {
    return entity.local_offset() != 0;
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
      const auto& hostPartition = hostPartitions[partitionIdx];
      const auto& partitionParts = hostPartition->get_legacy_partition_id();
      auto numBucketsInPartition = hostPartition->num_buckets();
      partitionBuffer.emplace_back(m_mesh, this, partitionParts, rank, partitionIdx);
      partitionBuffer[partitionIdx].m_buckets.reserve(numBucketsInPartition);

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
          bucketBuffer[iBucket].m_owningMesh = m_mesh;
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

        Kokkos::Profiling::pushRegion("partitionBuffer.add_bucket");
        // Update device buckets in partition
        partitionBuffer[partitionIdx].add_bucket(&bucketBuffer[iBucket]);
        Kokkos::Profiling::popRegion();

        STK_ThrowRequire(i == partitionBuffer[partitionIdx].num_buckets()-1);
        STK_ThrowRequire(iBucket == stkBucket.bucket_id());
      }
    }

    if (is_last_bucket_reference(rank)) {
      Kokkos::Profiling::pushRegion("dtors DevBkt");
      for (unsigned i = 0; i < m_buckets[rank].size(); ++i) {
        m_buckets[rank][i].~DeviceBucketT();
      }
      Kokkos::Profiling::popRegion();
    }

    if (is_last_partition_reference(rank)) {
      Kokkos::Profiling::pushRegion("dtors DevParttn");
      for (unsigned i = 0; i < m_partitions[rank].size(); ++i) {
        m_partitions[rank][i].~DevicePartition();
      }
      Kokkos::Profiling::popRegion();
    }

    m_buckets[rank] = bucketBuffer;
    m_partitions[rank] = partitionBuffer;

    copy_device_part_rank_info_from_host();
    Kokkos::Profiling::popRegion();
  }

  template<typename EntityViewType>
  unsigned get_max_num_parts_per_entity(EntityViewType const& entities) const;

  void sync_to_host(BucketRepository& hostBucketRepo);

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

  void check_mesh_consistency() const
  {
    STK_ThrowAssert(num_buckets(stk::topology::NODE_RANK) == m_buckets[stk::topology::NODE_RANK].size());
    STK_ThrowAssert(num_buckets(stk::topology::EDGE_RANK) == m_buckets[stk::topology::EDGE_RANK].size());
    STK_ThrowAssert(num_buckets(stk::topology::FACE_RANK) == m_buckets[stk::topology::FACE_RANK].size());
    STK_ThrowAssert(num_buckets(stk::topology::ELEM_RANK) == m_buckets[stk::topology::ELEM_RANK].size());

    STK_ThrowAssert(num_partitions(stk::topology::NODE_RANK) == m_partitions[stk::topology::NODE_RANK].size());
    STK_ThrowAssert(num_partitions(stk::topology::EDGE_RANK) == m_partitions[stk::topology::EDGE_RANK].size());
    STK_ThrowAssert(num_partitions(stk::topology::FACE_RANK) == m_partitions[stk::topology::FACE_RANK].size());
    STK_ThrowAssert(num_partitions(stk::topology::ELEM_RANK) == m_partitions[stk::topology::ELEM_RANK].size());
  }

  void print_part_ordinals() {
#ifndef NDEBUG
    std::cout << "DEBUGGING PRINTOUTS" << std::endl;
    for (auto rank = topology::BEGIN_RANK; rank < endRank; ++rank) {
      std::cout << num_partitions(rank) << " partitions for " << rank << std::endl;
      std::cout << num_buckets(rank) << " buckets for " << rank << std::endl;
    }

    std::cout << std::endl << "Printing Part Ordinals" << std::endl;
    for (auto rank = topology::BEGIN_RANK; rank < endRank; ++rank) {
      for (unsigned pid = 0; pid < m_partitions[rank].size(); ++pid) {
        auto& partition = m_partitions[rank][pid];
        if (!partition.is_active()) {
          std::cout << "\tPartition " << partition.partition_id() << "(stale) is ianctive" << std::endl;
          continue;
        }

        std::cout << "\tPartition [" << rank << "][" << partition.partition_id() << "] Part Ordinals: " << std::endl;
        auto hostPartOrdinals = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, partition.superset_part_ordinals());
        for (unsigned i = 0; i < hostPartOrdinals.extent(0); ++i) {
          std::cout << "\t" << hostPartOrdinals(i);
        }
        std::cout << std::endl;

        std::cout << "\t" << partition.num_buckets() << " buckets in this partition" << std::endl;

        for (unsigned i = 0; i < partition.m_buckets.size(); ++i) {
          auto bucket = partition.m_buckets[i];

          if (!bucket->is_active()) {
            std::cout << "\t\tBucket " << bucket->bucket_id() << "(stale) is inactive" << std::endl;
            continue;
          }
          std::cout << "\t\tBucket [" << rank << "][" << bucket->bucket_id() << "] Part Ordinals: " << std::endl;
          auto hostBucketPartOrdinals = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, bucket->get_part_ordinals());
          std::cout << "\t\t";
          for (unsigned j = 0; j < hostBucketPartOrdinals.extent(0); ++j) {
            std::cout << hostBucketPartOrdinals(j) << "\t";
          }
          std::cout << std::endl;

          std::cout << "\t\tThis bucket conatins " << bucket->size() << " entities:" << std::endl;
          auto hostEntities = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, bucket->m_entities);
          std::cout << "\t\t";
          for (unsigned k = 0; k < bucket->m_entities.extent(0); ++k) {
            std::cout << hostEntities(k) << " ";
          }
          std::cout << std::endl;
        }
      }
    }
#endif
  }

  unsigned m_initialBucketCapacity;
  unsigned m_maximumBucketCapacity;

  DeviceMeshT<NgpMemSpace>* m_mesh;

  DeviceBucketViewVector m_buckets[stk::topology::NUM_RANKS];
  DevicePartitionViewVector m_partitions[stk::topology::NUM_RANKS];

  // TODO Refactor to use a proper device side part (either DevicePart or a device copyable Part)
  DevicePartRankView m_partRanks;
  EntityRank endRank;
  unsigned m_numInternalParts;
  unsigned m_lastInternalPartOrdinal;

private:

  using SizeAndCapacity = DeviceFieldDataManagerBase::SizeAndCapacity;

  std::vector<unsigned> map_host_partitions_to_device_partitions(const BucketRepository& hostBucketRepo, EntityRank rank);

  void move_or_delete_partitions_on_host(BucketRepository& hostBucketRepo);

  void create_new_partitions_on_host(BucketRepository& hostBucketRepo);

  void create_new_buckets_on_host(BucketRepository& hostBucketRepo);

  void copy_partitions_and_buckets_to_host(BucketRepository& hostBucketRepo);

  void update_field_buckets(EntityRank rank, const std::vector<BucketShift>& bucketShifts);
};

template <typename NgpMemSpace>
template <typename EntityViewType, typename PartOrdinalProxyIndicesViewType, typename EntitySrcDestView>
void DeviceBucketRepository<NgpMemSpace>::batch_get_partitions(EntityViewType const& entities,
    PartOrdinalProxyIndicesViewType const& partOrdinalProxy,
    EntitySrcDestView& entitySrcDestView) const
{
  auto& deviceMesh = *m_mesh;
  Kokkos::parallel_for(
      entities.extent(0), KOKKOS_CLASS_LAMBDA(const int idx) {
        auto entity = entities(idx);
        auto rank = deviceMesh.entity_rank(entity);
        auto fastMeshIndex = deviceMesh.device_mesh_index(entity);
        auto srcBucketId = fastMeshIndex.bucket_id;
        auto srcBucketOrd = fastMeshIndex.bucket_ord;
        auto& bucket = deviceMesh.get_bucket(rank, srcBucketId);
        auto srcPartitionId = bucket.m_owningPartitionId;
        auto& partitions = m_partitions[rank];

        auto span = m_partitions[rank].size();

        entitySrcDestView(idx) = EntitySrcDest{
            entity, rank, srcPartitionId, srcBucketId, srcBucketOrd,
            INVALID_PARTITION_ID, INVALID_BUCKET_ID, INVALID_INDEX, INVALID_INDEX,
            0, 0, 0, 0};

        // TODO convert to a proper search after tuning into a nested team-based loop
        for (unsigned i = 0; i < span; ++i) {
          auto& partition = partitions[i];
          auto partitionId = partition.partition_id();
          if (partitionId == INVALID_PARTITION_ID) {
            continue;
          }

          auto partitionPartOrdinals = partition.superset_part_ordinals();
          if (partOrdinalProxy(idx).length != partitionPartOrdinals.extent(0)) {
            continue;
          }

          auto foundId = partitionId;
          for (unsigned j = 0; j < partOrdinalProxy(idx).length; j++) {
            auto localPartOrdinal = partOrdinalProxy(idx).startPtr + j;

            if (*localPartOrdinal != partitionPartOrdinals(j)) {
              foundId = INVALID_PARTITION_ID;
              break;
            }
          }
          if (foundId != INVALID_PARTITION_ID) {
            entitySrcDestView(idx).destPartitionId = foundId;
            break;
          }
        }
      });

  Kokkos::fence();
}

struct IsInvalidEntity
{
  KOKKOS_DEFAULTED_FUNCTION
  IsInvalidEntity() = default;

  KOKKOS_INLINE_FUNCTION
  bool operator()(Entity e) const {
    return e == Entity{};
  }
};

template <typename ConnectivityType>
struct IsInvalidConn
{
  KOKKOS_DEFAULTED_FUNCTION
  IsInvalidConn() = default;

  KOKKOS_INLINE_FUNCTION
  bool operator()(ConnectivityType c) const {
    return c == INVALID_CONNECTIVITY_ORDINAL;
  }
};

struct IsInvalidPerm
{
  KOKKOS_DEFAULTED_FUNCTION
  IsInvalidPerm() = default;

  KOKKOS_INLINE_FUNCTION
  bool operator()(Permutation p) const {
    return p == INVALID_PERMUTATION;
  }
};

template <typename NgpMemSpace>
void DeviceBucketRepository<NgpMemSpace>::compress_connectivities_in_buckets(EntityRank rank)
{
  using ExecSpace = typename NgpMemSpace::execution_space;
  using TeamHandle = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

  auto& buckets = m_buckets[rank];

  Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace>(buckets.size(), Kokkos::AUTO),
    KOKKOS_CLASS_LAMBDA(TeamHandle const& team) {
      auto bucketIdx = team.league_rank();
      auto& bucket = buckets[bucketIdx];

      if (!bucket.is_active() || bucket.size() == 0) { return; }

      // invalidate connectivities of invalid entities
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bucket.get_active_entity_span()),
        [&, this](const int bucketOrd) {
          auto entity = bucket[bucketOrd];

          if (is_valid(entity)) { return; }

          auto nodeConnOffset = bucket.get_node_connectivity_offset(bucketOrd);
          auto nodeConnLength = bucket.get_num_node_connectivity(bucketOrd);

          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, nodeConnOffset, nodeConnOffset+nodeConnLength),
            [&](const int connIdx) {
              bucket.m_nodeConnectivity(connIdx) = Entity{};
            }
          );

          for (auto lRank = stk::topology::BEGIN_RANK; lRank < endRank; lRank++) {
            auto sparseConnOffset = bucket.get_sparse_connectivity_offset(lRank, bucketOrd);
            auto sparseConnLength = bucket.get_num_sparse_connectivity(lRank, bucketOrd);

            Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, sparseConnOffset, sparseConnOffset+sparseConnLength),
              [&](const int connIdx) {
                bucket.m_sparseConnectivity(connIdx) = Entity{};
                bucket.m_sparseConnectivityOrdinals(connIdx) = INVALID_CONNECTIVITY_ORDINAL;
                if (bucket.has_permutations()) {
                  bucket.m_sparseConnectivityPermutations(connIdx) = INVALID_PERMUTATION;
                }
              }
            );
          }
        }
      );

      team.team_barrier();

      // temporarily replace offsets in offset views with each connectivity length
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bucket.get_active_entity_span()),
        [&, this](const int bucketOrd) {
          auto entity = bucket[bucketOrd];

          if (!is_valid(entity)) { return; }

          auto nodeConnLength = bucket.get_num_node_connectivity(bucketOrd);

          bucket.m_nodeConnectivityOffsets(bucketOrd).offset = nodeConnLength;

          for (auto lRank = stk::topology::EDGE_RANK; lRank < endRank; lRank++) {
            auto sparseConnLength = bucket.get_num_sparse_connectivity(lRank, bucketOrd);
            bucket.m_sparseConnectivityOffsets(lRank, bucketOrd).offset = sparseConnLength;
          }
        }
      );

      team.team_barrier();

      // remove connectivities with invalid ordinals
      Kokkos::Experimental::remove_if(team, bucket.m_nodeConnectivity, IsInvalidEntity{});
      Kokkos::Experimental::remove_if(team, bucket.m_sparseConnectivity, IsInvalidEntity{});
      Kokkos::Experimental::remove_if(team, bucket.m_sparseConnectivityOrdinals, IsInvalidConn<typename decltype(bucket.m_sparseConnectivityOrdinals)::value_type>{});
      Kokkos::Experimental::remove_if(team, bucket.m_sparseConnectivityPermutations, IsInvalidPerm{});

      team.team_barrier();

      // fix offsets based on "length"
      Kokkos::single(Kokkos::PerTeam(team),
        [&, this]() {
          auto nodeOffset = 0;
          auto sparseOffset = 0;
          for (unsigned i = 0; i < bucket.get_active_entity_span(); ++i) {
            if (bucket[i] == Entity{}) { continue; }

            // temporarily stored as length
            auto nodeLength = bucket.m_nodeConnectivityOffsets(i).offset;
            bucket.m_nodeConnectivityOffsets(i).offset = nodeOffset;
            nodeOffset += nodeLength;

            auto sparseLength = 0u;
            for (auto lRank = stk::topology::EDGE_RANK; lRank < endRank; lRank++) {
              sparseLength = bucket.m_sparseConnectivityOffsets(lRank, i).offset;
              bucket.m_sparseConnectivityOffsets(lRank, i).offset = sparseOffset;
              sparseOffset += sparseLength;
            }
          }

          unsigned lastOffset = bucket.get_active_entity_span();
          bucket.m_nodeConnectivityOffsets(lastOffset).offset = nodeOffset;
          bucket.m_sparseConnectivityOffsets(endRank+1, bucket.get_active_entity_span()-1).offset = sparseOffset;
        }
      );
    }
  );
  Kokkos::fence();
}

// TODO use Kokkos::copy(team) when batch_move_entities() is refactored
template <typename DeviceBucket>
KOKKOS_FUNCTION void copy_bucket_connectivities_from_src(
    DeviceBucket* srcBucket, DeviceBucket* destBucket, EntitySrcDest const& entitySrcDest, EntityRank endRank)
{
  auto srcBucketOrd = entitySrcDest.srcBucketOrd;
  auto destBucketOrd = entitySrcDest.destBucketOrd;

  auto destNodeConnectivityStartIdx = entitySrcDest.destNodeConnectivityStartIdx;
  auto destNodeConnectivityLength = entitySrcDest.destNodeConnectivityLength;

  auto destSparseConnectivityStartIdx = entitySrcDest.destSparseConnectivityStartIdx;
  auto destSparseConnectivityLength = entitySrcDest.destSparseConnectivityLength;

#ifndef NDEBUG
  auto srcNumNodeConnectivityLength = srcBucket->get_num_node_connectivity(srcBucketOrd);
  if (srcNumNodeConnectivityLength != Kokkos::Experimental::finite_max_v<ConnectivityOrdinal> &&
      srcNumNodeConnectivityLength != destNodeConnectivityLength) {
    Kokkos::abort("srcNumNodeConnectivityLength does not match destNodeConnectivityLength");
  }

  unsigned totalSrcNumSparseConnectivityLength = 0;
  for (auto rank = stk::topology::EDGE_RANK; rank < endRank; ++rank) {
    unsigned srcNumSparseConnectivityLength = srcBucket->get_num_sparse_connectivity(rank, srcBucketOrd);

    if (srcNumSparseConnectivityLength == Kokkos::Experimental::finite_max_v<ConnectivityOrdinal>) {
      continue;
    }

    totalSrcNumSparseConnectivityLength += srcNumSparseConnectivityLength;
  }
  if (totalSrcNumSparseConnectivityLength != destSparseConnectivityLength) {
    Kokkos::abort("srcNumSparseConnectivityLength does not match destSparseConnectivityLength");
  }
#endif

  // copy node connectivity offset
  destBucket->m_nodeConnectivityOffsets(destBucketOrd) = {
      entitySrcDest.entity, static_cast<ConnectivityOrdinal>(destNodeConnectivityStartIdx)};

  // FIXME placing in last offset to correctly compute length
  if (destBucketOrd+2 == destBucket->m_nodeConnectivityOffsets.extent(0)) {
    destBucket->m_nodeConnectivityOffsets(destBucketOrd+1) = {
        entitySrcDest.entity, static_cast<ConnectivityOrdinal>(destNodeConnectivityStartIdx+destNodeConnectivityLength)};
  }

  auto srcNodeConnectivities = srcBucket->get_connected_entities(srcBucketOrd, stk::topology::NODE_RANK);

  // copy node connectivity
  for (unsigned i = 0; i < destNodeConnectivityLength; ++i) {
    auto offset = destNodeConnectivityStartIdx + i;
    destBucket->m_nodeConnectivity(offset) = srcNodeConnectivities[i];
  }

  auto srcSparseConnectivityStartIdx = srcBucket->m_sparseConnectivityOffsets(topology::EDGE_RANK, srcBucketOrd).offset;
#ifndef NDEBUG
  if (srcSparseConnectivityStartIdx == INVALID_INDEX) {
    Kokkos::abort("srcSparseConnectivityOffset is invalid");
  }
#endif

  auto initOffset = destSparseConnectivityStartIdx;
  for (auto rank = stk::topology::EDGE_RANK; rank < endRank; ++rank) {
    auto length = srcBucket->get_num_sparse_connectivity(rank, srcBucketOrd);
    destBucket->m_sparseConnectivityOffsets(rank, destBucketOrd) = impl::EntityOffsetComp<unsigned>{entitySrcDest.entity, initOffset};
    initOffset += length;
  }
  destBucket->m_sparseConnectivityOffsets(endRank, destBucketOrd) = impl::EntityOffsetComp<unsigned>{entitySrcDest.entity, initOffset};

  // WITH continuous mmeory assumption
  for (unsigned i = 0; i < destSparseConnectivityLength; ++i) {
    auto destOffset = destSparseConnectivityStartIdx + i;
    auto srcSparseConnectivityOffset  = srcSparseConnectivityStartIdx + i;
    destBucket->m_sparseConnectivity(destOffset) = srcBucket->m_sparseConnectivity(srcSparseConnectivityOffset);
    destBucket->m_sparseConnectivityOrdinals(destOffset) = srcBucket->m_sparseConnectivityOrdinals(srcSparseConnectivityOffset);
#ifndef NDEBUG
    if (srcBucket->has_permutations() != destBucket->has_permutations()) {
      Kokkos::abort("src and dest bucket do not have matching permutation property");
    }
#endif
    if (srcBucket->has_permutations()) {
      destBucket->m_sparseConnectivityPermutations(destOffset) = srcBucket->m_sparseConnectivityPermutations(srcSparseConnectivityOffset);
    }
  }
}


template <typename NgpMemSpace>
void DeviceBucketRepository<NgpMemSpace>::sync_to_host(BucketRepository& hostBucketRepo)
{
  create_new_partitions_on_host(hostBucketRepo);
  create_new_buckets_on_host(hostBucketRepo);
  copy_partitions_and_buckets_to_host(hostBucketRepo);
}

template <typename NgpMemSpace>
std::vector<unsigned> DeviceBucketRepository<NgpMemSpace>::map_host_partitions_to_device_partitions(const BucketRepository& hostBucketRepo, EntityRank rank)
{
  DevicePartitionViewVector& devicePartitions = get_partitions(rank);
  std::vector<std::vector<stk::mesh::PartOrdinal>> devicePartOrdinals;
  for (unsigned i=0; i < devicePartitions.size(); ++i)
  {
    auto ordinals = devicePartitions[i].superset_part_ordinals();
    auto ordinalsHost = Kokkos::create_mirror_view(ordinals);
    Kokkos::deep_copy(ordinalsHost, ordinals);
    std::vector<stk::mesh::PartOrdinal>& ordinalsVec = devicePartOrdinals.emplace_back();
    for (unsigned j=0; j < ordinalsHost.extent(0); ++j)
    {
      ordinalsVec.push_back(ordinalsHost[j]);
    }
  }

  constexpr unsigned INVALID = std::numeric_limits<unsigned>::max();
  const std::vector<Partition*>& hostPartitions = hostBucketRepo.m_partitions[rank];
  std::vector<unsigned> hostPartitionIdxOnDevice(hostPartitions.size());
  for (unsigned i=0; i < hostPartitions.size(); ++i)
  {
    const std::vector<stk::mesh::PartOrdinal> hostPartOrdinals = hostPartitions[i]->get_legacy_partition_id();
    auto it = std::lower_bound(devicePartOrdinals.begin(), devicePartOrdinals.end(), hostPartOrdinals, DevicePartOrdinalLess{});
    if (*it == hostPartOrdinals)
    {
      hostPartitionIdxOnDevice[i] = std::distance(devicePartOrdinals.begin(), it);
    } else
    {
      hostPartitionIdxOnDevice[i] = INVALID;
    }
  }

  return hostPartitionIdxOnDevice;
}

template <typename NgpMemSpace>
void DeviceBucketRepository<NgpMemSpace>::move_or_delete_partitions_on_host(BucketRepository& hostBucketRepo)
{
  constexpr unsigned INVALID = std::numeric_limits<unsigned>::max();
  for (EntityRank rank=stk::topology::NODE_RANK; rank < hostBucketRepo.mesh().mesh_meta_data().entity_rank_count(); ++rank)
  {
    std::vector<unsigned> hostPartitionIdxOnDevice = map_host_partitions_to_device_partitions(hostBucketRepo, rank);
    std::vector<Partition*> permutedPartitions(num_partitions(rank), nullptr);
    std::vector<Partition*>& hostPartitions = hostBucketRepo.m_partitions[rank];

    for (unsigned i=0; i < hostPartitionIdxOnDevice.size(); ++i)
    {
      if (hostPartitionIdxOnDevice[i] != INVALID)
      {
        permutedPartitions[hostPartitionIdxOnDevice[i]] = hostPartitions[i];
      } else
      {
        hostBucketRepo.deallocate_partition(rank, i);
      }
    }

    hostPartitions = permutedPartitions;
  }
}


template <typename NgpMemSpace>
void DeviceBucketRepository<NgpMemSpace>::create_new_partitions_on_host(BucketRepository& hostBucketRepo)
{
  move_or_delete_partitions_on_host(hostBucketRepo);

  for (EntityRank rank=stk::topology::NODE_RANK; rank < hostBucketRepo.mesh().mesh_meta_data().entity_rank_count(); ++rank)
  {
    STK_ThrowAssert(num_partitions(rank) == hostBucketRepo.get_partitions(rank).size());

    DevicePartitionViewVector& devicePartitions = get_partitions(rank);
    for (unsigned i=0; i < num_partitions(rank); ++i)
    {
      if (!hostBucketRepo.get_partition(rank, i))
      {
        DevicePartition<NgpMemSpace>& devicePartition = devicePartitions[i];

        const PartOrdinalViewType<NgpMemSpace>& deviceOrdinals = devicePartition.superset_part_ordinals();
        auto hostOrdinals = Kokkos::create_mirror_view(deviceOrdinals);
        Kokkos::deep_copy(hostOrdinals, deviceOrdinals);
        OrdinalVector ordinalsVec(hostOrdinals.extent(0));
        for (unsigned j=0; j < hostOrdinals.extent(0); ++j)
        {
          ordinalsVec[j] = hostOrdinals(j);
        }

        auto i_iterator = hostBucketRepo.m_partitions[rank].begin() + i;
        hostBucketRepo.create_partition(rank, ordinalsVec, i_iterator);
      }
    }
  }
}

template <typename NgpMemSpace>
void DeviceBucketRepository<NgpMemSpace>::create_new_buckets_on_host(BucketRepository& hostBucketRepo)
{
  // create empty buckets on host, such that Partition::get_bucket(i) and DevicePartitions::get_bucket(i) give
  // corresponding buckets
  for (EntityRank rank=stk::topology::NODE_RANK; rank < hostBucketRepo.mesh().mesh_meta_data().entity_rank_count(); ++rank)
  {
    const std::vector<Partition*> partitions = hostBucketRepo.get_partitions(rank);
    DevicePartitionViewVector& devicePartitions = get_partitions(rank);
    STK_ThrowRequireMsg(partitions.size() == devicePartitions.size(),
                       "number of host and device partitions should match");
    for (unsigned p=0; p < devicePartitions.size(); ++p)
    {
      Partition& hostPartition                       = *(partitions[p]);
      DevicePartition<NgpMemSpace> & devicePartition = devicePartitions[p];

      if (hostPartition.num_buckets() > devicePartition.num_buckets())
      {
        STK_ThrowRequireMsg(devicePartition.num_buckets() > 0, "empty partitions should have been removed during meshmod");
        for (unsigned i=hostPartition.num_buckets()-1; i >= devicePartition.num_buckets(); --i)
        {
          hostPartition.delete_bucket(hostPartition.get_bucket(i));
        }
      } else
      {
        for (unsigned i=hostPartition.num_buckets(); i < devicePartition.num_buckets(); ++i)
        {
          hostPartition.add_empty_bucket();
        }
      }
    }
  }
}

template <typename NgpMemSpace>
void DeviceBucketRepository<NgpMemSpace>::copy_partitions_and_buckets_to_host(BucketRepository& hostBucketRepo)
{
  for (EntityRank rank=stk::topology::NODE_RANK; rank < hostBucketRepo.mesh().mesh_meta_data().entity_rank_count(); ++rank)
  {
    const DevicePartitionViewVector& devicePartitions = get_partitions(rank);
    for (unsigned i=0; i < devicePartitions.size(); ++i)
    {
      Partition* hostPartition = hostBucketRepo.get_partition(rank, i);
      DevicePartition<NgpMemSpace>& devicePartition = devicePartitions[i];
      devicePartition.sync_to_host(*hostPartition);
    }
  }
}

template <typename NgpMemSpace>
void DeviceBucketRepository<NgpMemSpace>::update_field_buckets(EntityRank rank, const std::vector<BucketShift>& bucketShifts)
{
  std::vector<SizeAndCapacity> sizes;

  sizes.clear();
  for (unsigned i=0; i < num_buckets(rank); ++i)
  {
    const DeviceBucket& bucket = *(get_bucket(rank, i));
    sizes.push_back(SizeAndCapacity{bucket.size(), bucket.capacity()});
  }

  DeviceFieldDataManagerBase* fieldDataManager = m_mesh->get_field_data_manager(m_mesh->get_bulk_on_host());

  fieldDataManager->reorder_and_resize_buckets(rank, m_mesh->get_bulk_on_host().mesh_meta_data().get_fields(rank), sizes, bucketShifts);
}

}  // namespace impl
}  // namespace mesh
}  // namespace stk

#endif
