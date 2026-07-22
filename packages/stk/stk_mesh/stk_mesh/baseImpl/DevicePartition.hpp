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

#ifndef STK_MESH_DEVICE_PARTITION_HPP
#define STK_MESH_DEVICE_PARTITION_HPP

#include "Kokkos_Atomic.hpp"
#include "Kokkos_Macros.hpp"
#include "Kokkos_Core.hpp"
#include "stk_util/stk_config.h"
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/base/DeviceBucket.hpp"
#include "stk_mesh/baseImpl/DeviceMeshViewVector.hpp"
#include "stk_mesh/baseImpl/NgpMeshImpl.hpp"
#include "stk_mesh/baseImpl/ViewVector.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

namespace stk {
namespace mesh {

template<typename NgpMemSpace>
class DeviceMeshT;

namespace impl {

template <typename NgpMemSpace>
class DeviceBucketRepository;

template <typename NgpMemSpace>
struct DeviceBucketPtrWrapper
{
  using Type = DeviceBucketT<NgpMemSpace>;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr DeviceBucketPtrWrapper() = default;

  KOKKOS_FUNCTION
  DeviceBucketPtrWrapper(Type* ptr) : bucketPtr(ptr), bucketId(ptr->bucket_id()) {}

  KOKKOS_FUNCTION
  Type& operator*() {
    return *bucketPtr;
  }

  KOKKOS_FUNCTION
  const Type& operator*() const {
    return *bucketPtr;
  }
  
  KOKKOS_FUNCTION
  Type* operator->() const {
    return bucketPtr;
  }

  KOKKOS_FUNCTION
  bool operator<(DeviceBucketPtrWrapper<NgpMemSpace> const& other) const {
    return *bucketPtr < *other.bucketPtr;
  }

  Type* bucketPtr = nullptr;
  unsigned bucketId = INVALID_BUCKET_ID;
};

template <typename NgpMemSpace>
class DevicePartition
{
 public:
  using DeviceBucket = DeviceBucketT<NgpMemSpace>;
  using DeviceBucketPtrViewVector = ImplDeviceMeshViewVector<DeviceBucketPtrWrapper<NgpMemSpace>, stk::ngp::UVMMemSpace>;
  using DeviceBucketPtrUView = Kokkos::View<DeviceBucketPtrWrapper<NgpMemSpace>*, stk::ngp::UVMMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  KOKKOS_FUNCTION
  DevicePartition()
    : m_mesh(nullptr),
      m_deviceBucketRepo(nullptr),
      m_rank(stk::topology::INVALID_RANK),
      m_partitionId(INVALID_PARTITION_ID),
      m_numEntities(0),
      m_isModified(false)
  {}

  DevicePartition(DeviceMeshT<NgpMemSpace>* deviceMesh, DeviceBucketRepository<NgpMemSpace>* deviceBucketRepo,
                  PartOrdinalViewType<NgpMemSpace> const& partOrdinals, EntityRank rank, unsigned partitionId,
                  bool copyPartOrdinals = false)
    : m_mesh(deviceMesh),
      m_deviceBucketRepo(deviceBucketRepo),
      m_partOrdinals(partOrdinals),
      m_buckets("DeviceBuckets"),
      m_rank(rank),
      m_partitionId(partitionId),
      m_numEntities(0),
      m_isModified(false)
  {
    if (copyPartOrdinals) {
      m_partOrdinals = PartOrdinalViewType<NgpMemSpace>("DevicePartitionPartOrdinals", partOrdinals.extent(0));
      Kokkos::deep_copy(m_partOrdinals, partOrdinals);
    }
  }

  DevicePartition(DeviceMeshT<NgpMemSpace>* deviceMesh, DeviceBucketRepository<NgpMemSpace>* deviceBucketRepo,
                  std::vector<PartOrdinal> const& partOrdinals, EntityRank rank, unsigned partitionId)
    : m_mesh(deviceMesh),
      m_deviceBucketRepo(deviceBucketRepo),
      m_buckets("DeviceBuckets"),
      m_rank(rank),
      m_partitionId(partitionId),
      m_numEntities(0),
      m_isModified(false)
  {
    set_part_ordinals_from_host(partOrdinals);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned partition_id() const { return m_partitionId; }

  void set_partition_id(unsigned id) { m_partitionId = id; }

  KOKKOS_INLINE_FUNCTION
  EntityRank get_rank() const { return m_rank; }

  KOKKOS_INLINE_FUNCTION
  bool is_active() const { return m_partitionId != INVALID_PARTITION_ID; }

  KOKKOS_INLINE_FUNCTION
  bool is_modified() const { return m_isModified; }

  KOKKOS_INLINE_FUNCTION
  void set_modified() { m_isModified = true; }

  KOKKOS_INLINE_FUNCTION
  void clear_modified() { m_isModified = false; }

  void sync_to_host(Partition& hostPartition) const
  {
#ifndef NDEBUG
    STK_ThrowRequireMsg(get_rank() == hostPartition.get_rank(),
                        "Cannot copy DeviceParition to HostParition of different entity rank");
    auto partOrdinalsFromDevice = Kokkos::create_mirror_view(m_partOrdinals);
    Kokkos::deep_copy(partOrdinalsFromDevice, m_partOrdinals);
    STK_ThrowRequireMsg(partOrdinalsFromDevice.extent(0), hostPartition.get_legacy_partition_id().size());
    for (size_t i=0; i < partOrdinalsFromDevice.extent(0); ++i)
    {
      STK_ThrowRequireMsg(partOrdinalsFromDevice(i) == hostPartition.get_legacy_partition_id()[i],
                          "Cannot copy DevicePartition to HostPartition with different part ordinals");
    }
#endif

    for (size_t i=0; i < num_buckets(); ++i)
    {
      DeviceBucketPtrWrapper<NgpMemSpace>& deviceBucket = m_buckets[i];
      Bucket* hostBucket = hostPartition.get_bucket(i);
      deviceBucket->sync_to_host(*hostBucket);
    }
    hostPartition.compute_size();
  }

  KOKKOS_INLINE_FUNCTION
  bool add_entity_without_connectivity(const Entity entity, unsigned destBucketIndex, unsigned destBucketOrd)
  {
    if (destBucketIndex >= m_buckets.size()) {
      Kokkos::abort("Bucket index larger than bucket size");
    }
    m_buckets[destBucketIndex]->add_entity(entity, destBucketOrd);

    update_partition_meta_entity_added();
    return true;
  }

  KOKKOS_FUNCTION
  bool remove_entity(unsigned srcBucketId, unsigned srcBucketOrd)
  {
    auto bucketOrdInPartition = INVALID_INDEX;
    for (unsigned i = 0; i < num_buckets(); ++i) {
      if (m_buckets[i].bucketId == srcBucketId) {
        bucketOrdInPartition = i;
        break;
      }
    }

    if (bucketOrdInPartition == INVALID_INDEX) {
      Kokkos::abort("Matching bucket not found");
    }
    m_buckets[bucketOrdInPartition]->remove_entity(srcBucketOrd);

    update_partition_meta_entity_removed();
    return true;
  }

  void add_bucket(DeviceBucket* bucket)
  {
    bucket->m_owningPartitionId = m_partitionId;

    m_buckets.emplace_back(bucket);

    update_partition_meta_bucket_added(bucket);
  }

  void remove_bucket(DeviceBucket* bucket)
  {
    update_partition_meta_bucket_removed(bucket);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_last_avail_bucket_index() const
  {
    if (num_buckets() == 0) {
      return INVALID_INDEX;
    }

    unsigned lastBucketIdx = m_buckets.size() - 1;
    return (!m_buckets[lastBucketIdx]->is_active() || m_buckets[lastBucketIdx]->is_full()) ? INVALID_INDEX
                                                                                           : lastBucketIdx;
  }

  KOKKOS_INLINE_FUNCTION
  DeviceBucketPtrWrapper<NgpMemSpace>& get_bucket(unsigned idx) const
  {
    return m_buckets[idx];
  }

  void sort_buckets()
  {
    using ExecSpace = typename NgpMemSpace::execution_space;
    DeviceBucketPtrUView compactBucketPtrUView(m_buckets.data(), m_buckets.size());
    if (!Kokkos::Experimental::is_sorted(ExecSpace{}, compactBucketPtrUView)) {
      Kokkos::sort(ExecSpace{}, compactBucketPtrUView);
    }
    STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, compactBucketPtrUView));
    m_buckets.resize(num_buckets());
  }

  template <typename TeamMember, typename EntityView>
  KOKKOS_INLINE_FUNCTION
  void gather_all_entities_and_sort(TeamMember const& teamMember, EntityView& allEntities)
  {
    auto bucketCapacity = get_bucket(0)->capacity();

    // gather all entities
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_buckets()),
      [&](const int bucketIdx) {
        auto bucket = get_bucket(bucketIdx);
        auto isActive = bucket->is_active() && bucket->size() > 0;
        auto destIdx = bucketIdx * bucketCapacity;
        auto& srcEntities = bucket->m_entities;
        auto bucketActiveSpan = bucket->get_active_entity_span();

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, bucketCapacity),
          [&](const int entityIdx) {
            auto invalid = !isActive || (static_cast<unsigned>(entityIdx) >= bucketActiveSpan);
            allEntities(destIdx+entityIdx) = invalid ? Entity{} : srcEntities(entityIdx);
          }
        );
      }
    );
    teamMember.team_barrier();

    if (!Kokkos::Experimental::is_sorted(teamMember, allEntities, EntityCompareInvalidAtEnd{})) {
      Kokkos::Experimental::sort_team(teamMember, allEntities, EntityCompareInvalidAtEnd{});
    }

#ifndef NDEBUG
    auto sorted = Kokkos::Experimental::is_sorted(teamMember, allEntities, EntityCompareInvalidAtEnd{});

    // check that gathered valid entity count equals entity count known by this partition
    Kokkos::single(Kokkos::PerTeam(teamMember), [&]() {
      STK_NGP_ThrowRequireMsg(sorted, "Entities are not sorted");

      unsigned countEntities = 0;
      for (unsigned i = 0; i < allEntities.extent(0); ++i) {
        if (allEntities[i].m_value != Entity::InvalidEntity) countEntities++;
      }
      STK_NGP_ThrowRequireMsg(countEntities == num_entities(), "Number of entities in a partition does not match active entities in its buckets.");
    });
#endif
    teamMember.team_barrier();
  }

  template <typename TeamMember, typename EntityView>
  KOKKOS_INLINE_FUNCTION
  void scatter_entities_to_buckets(TeamMember const& teamMember, EntityView& allEntities)
  {
    auto bucketCapacity = get_bucket(0)->capacity();

    // scatter to buckets and update bucket info
    auto numBucketsNeeded = (num_entities() + bucketCapacity - 1) / bucketCapacity;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_buckets()),
      [&](const unsigned bucketIdx) {
        auto bucket = get_bucket(bucketIdx);

        if (bucketIdx >= numBucketsNeeded) {
          bucket->clear_entities();
          return;
        }

        auto srcIdx = bucketIdx * bucket->capacity();
        auto& destEntities = bucket->m_entities;

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, bucketCapacity), 
          [&](const int entityIdx) {
            destEntities(entityIdx) = allEntities(srcIdx+entityIdx);
          }
        );

        auto getEntityCount = [&](int idx) { return Kokkos::min(bucketCapacity, num_entities() - idx * bucketCapacity); };

        Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
          bucket->m_bucketSize = getEntityCount(bucketIdx);
          bucket->m_activeEntitySpan = bucket->m_bucketSize;
          bucket->set_modified();
        });
      }
    );
    teamMember.team_barrier();

    // update counts in this partition
    Kokkos::single(Kokkos::PerTeam(teamMember), [&]() {
      m_buckets.set_active_entries(numBucketsNeeded);
      set_modified();
    });
  }

  template <typename TeamMember>
  KOKKOS_FUNCTION
  void sort_entities_in_buckets(TeamMember const& teamMember, unsigned teamScratchLevel, size_t perTeamScratchSize)
  {
    using ScratchView = Kokkos::View<Entity*, typename TeamMember::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    auto entityViewExtent = perTeamScratchSize / sizeof(Entity);
    ScratchView allEntities(teamMember.team_scratch(teamScratchLevel), entityViewExtent);

    gather_all_entities_and_sort(teamMember, allEntities);

    scatter_entities_to_buckets(teamMember, allEntities);
  }

  template <typename TeamMember, typename EntityView>
  requires (Kokkos::is_view_v<EntityView>)
  KOKKOS_FUNCTION
  void sort_entities_in_buckets(TeamMember const& teamMember, EntityView& allEntities, unsigned maxNumEntityPerPartition)
  {
    using EntityUView = Kokkos::View<Entity*, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    auto partitionIdx = teamMember.league_rank();
    auto startIdx = partitionIdx * maxNumEntityPerPartition;

    EntityUView entitiesInPartition(allEntities.data()+startIdx, maxNumEntityPerPartition);

    gather_all_entities_and_sort(teamMember, entitiesInPartition);

    scatter_entities_to_buckets(teamMember, entitiesInPartition);
  }

  void set_part_ordinals_from_host(std::vector<PartOrdinal> const& partOrdinals)
  {
    m_partOrdinals = PartOrdinalViewType<NgpMemSpace>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PartOrdinals"), partOrdinals.size());
    auto hostPartOrdinals = HostPartOrdinalViewType(partOrdinals.data(), partOrdinals.size());
    Kokkos::deep_copy(m_partOrdinals, hostPartOrdinals);
  }

  void reset_partition_id_in_owned_buckets()
  {
    for (unsigned i = 0; i < num_buckets(); ++i) {
      auto bucket = m_buckets[i].bucketPtr;
      bucket->m_owningPartitionId = m_partitionId;
    }
  }

  void update_partition_meta_bucket_added(DeviceBucket* addedBucket)
  {
    m_numEntities += addedBucket->size();
    set_modified();
  }

  void update_partition_meta_bucket_removed(DeviceBucket* removedBucket)
  {
    m_buckets.decrement_num_active_entries();
    m_numEntities -= removedBucket->size();
    removedBucket->m_owningPartitionId = INVALID_PARTITION_ID;
    set_modified();
  }

  KOKKOS_FUNCTION
  void update_partition_meta_entity_added()
  {
    Kokkos::atomic_inc(&m_numEntities);
    set_modified();
  }

  KOKKOS_FUNCTION
  void update_partition_meta_entity_removed()
  {
    Kokkos::atomic_dec(&m_numEntities);
    set_modified();
  }

  KOKKOS_INLINE_FUNCTION
  const PartOrdinalViewType<NgpMemSpace>& superset_part_ordinals() const { return m_partOrdinals; }

  KOKKOS_INLINE_FUNCTION
  size_t num_buckets() const { return m_buckets.num_active_entries(); }

  KOKKOS_INLINE_FUNCTION
  size_t num_entities() const { return m_numEntities; }

  KOKKOS_INLINE_FUNCTION
  bool has_no_buckets() const { return m_buckets.num_active_entries() == 0; }

  KOKKOS_INLINE_FUNCTION
  bool operator<(DevicePartition<NgpMemSpace> const& other) const {
    if (m_partitionId == INVALID_PARTITION_ID) {
      return false;
    } else if (other.m_partitionId == INVALID_PARTITION_ID) {
      return true;
    } else {
      return DevicePartOrdinalLess{}(superset_part_ordinals(), other.superset_part_ordinals());
    }
  }

  KOKKOS_FUNCTION
  void update_bucket_ptrs(DeviceBucketT<NgpMemSpace>& newBucket)
  {
    auto newBucketId = newBucket.bucket_id();
    auto bucketSpan = m_buckets.size();

    for (unsigned i = 0; i < bucketSpan; ++i) {
      auto oldBucketId = m_buckets[i].bucketId;

      if (oldBucketId == INVALID_BUCKET_ID || newBucketId == INVALID_BUCKET_ID) {
        continue;
      }

      if (oldBucketId == newBucketId) {
        m_buckets[i].bucketPtr = &newBucket;
        return;
      } else {
        // newly constructed bucket
      }
    }
  }
  
  DeviceMeshT<NgpMemSpace>* m_mesh;
  DeviceBucketRepository<NgpMemSpace>* m_deviceBucketRepo;
  PartOrdinalViewType<NgpMemSpace> m_partOrdinals;
  DeviceBucketPtrViewVector m_buckets;
  EntityRank m_rank;
  unsigned m_partitionId;
  unsigned m_numEntities;
  bool m_isModified;
};

template <typename SearchView, typename PartOrdinalView>
struct DeviceBucketPartsMatcherFunctor {
  using value_type = unsigned;

  KOKKOS_INLINE_FUNCTION
  DeviceBucketPartsMatcherFunctor(SearchView const& view, PartOrdinalView const& partOrdinals)
    : m_searchView(view),
      m_partOrdinalView(partOrdinals)
  {}

  KOKKOS_INLINE_FUNCTION
  void init(value_type& init) const {
    init = INVALID_PARTITION_ID;
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, value_type const& src) const {
    if (dst != src) {
      dst = (src == INVALID_PARTITION_ID) ? dst : src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void final(value_type& /*dst*/) const {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type& update) const {
    if (m_partOrdinalView.extent(0) == 0) { return; }
    if (m_searchView[i].partition_id() == INVALID_PARTITION_ID) { return; }

    auto partOrdinals = m_searchView[i].superset_part_ordinals();

    if (m_partOrdinalView.extent(0) != partOrdinals.extent(0)) { return; }

    for (unsigned j = 0; j < m_partOrdinalView.extent(0); ++j) {
      if (partOrdinals(j) != m_partOrdinalView(j)) {
        return;
      }
    }
    update = m_searchView[i].partition_id();
  }

  SearchView m_searchView;
  PartOrdinalView m_partOrdinalView;
};

template <typename PolicyType, typename PartitionView, typename PartOrdinalView>
unsigned
search_matching_device_partitions(PolicyType const& policy, PartitionView const& partitions, PartOrdinalView const& partOrdinals)
{
  DeviceBucketPartsMatcherFunctor functor(partitions, partOrdinals);
  unsigned searchResult = INVALID_PARTITION_ID;

  Kokkos::parallel_reduce(policy, functor, searchResult);
  return searchResult;
}

} } }

#endif
