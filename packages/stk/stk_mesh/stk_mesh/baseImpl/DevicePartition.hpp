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
#include "stk_mesh/baseImpl/NgpMeshImpl.hpp"
#include "Kokkos_Core.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace mesh {

template<typename NgpMemSpace>
class DeviceMeshT;

namespace impl {

constexpr unsigned initialDeviceBucketViewCapacity = 32;
constexpr unsigned initialDevicePartitionViewCapacity = 32;
constexpr unsigned INVALID_INDEX = std::numeric_limits<unsigned>::max();

template <typename NgpMemSpace>
class DeviceBucketRepository;

template <typename NgpMemSpace>
struct DeviceBucketWrapper
{
  using Type = DeviceBucketT<NgpMemSpace>;

  DeviceBucketWrapper() = default;

  DeviceBucketWrapper(Type* ptr)
    : bucketPtr(ptr),
      m_localBucketId(bucketPtr->bucket_id())
  {}

  operator Type() const {
    return *bucketPtr;
  }
  
  KOKKOS_FUNCTION
  bool operator<(DeviceBucketWrapper<NgpMemSpace> const& other) const {
    return m_localBucketId < other.m_localBucketId;
  }
  
  Type* bucketPtr = nullptr;
  unsigned m_localBucketId = INVALID_BUCKET_ID;
};

template <typename NgpMemSpace>
class DevicePartition
{
 public:
  using DeviceBucket = DeviceBucketT<NgpMemSpace>;
  using DeviceBucketPtrView = Kokkos::View<DeviceBucketWrapper<NgpMemSpace>*, stk::ngp::UVMMemSpace>;
  using DeviceBucketPtrUView = Kokkos::View<DeviceBucketWrapper<NgpMemSpace>*, stk::ngp::UVMMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  KOKKOS_FUNCTION
  DevicePartition()
    : m_mesh(nullptr),
      m_deviceBucketRepo(nullptr),
      m_rank(stk::topology::NODE_RANK),
      m_partitionId(INVALID_PARTITION_ID),
      m_numBuckets(0),
      m_bucketViewLastBucketIdx(INVALID_INDEX),
      m_numEntities(0)
  {}

  DevicePartition(DeviceMeshT<NgpMemSpace>* deviceMesh, DeviceBucketRepository<NgpMemSpace>* deviceBucketRepo,
                  PartOrdinalViewType<NgpMemSpace> const& partOrdinals, EntityRank rank, unsigned partitionId)
    : m_mesh(deviceMesh),
      m_deviceBucketRepo(deviceBucketRepo),
      m_partOrdinals(partOrdinals),
      m_rank(rank),
      m_partitionId(partitionId),
      m_numBuckets(0),
      m_bucketViewLastBucketIdx(INVALID_INDEX),
      m_numEntities(0)
  {}

  DevicePartition(DeviceMeshT<NgpMemSpace>* deviceMesh, DeviceBucketRepository<NgpMemSpace>* deviceBucketRepo,
                  std::vector<PartOrdinal> const& partOrdinals, EntityRank rank, unsigned partitionId)
    : m_mesh(deviceMesh),
      m_deviceBucketRepo(deviceBucketRepo),
      m_rank(rank),
      m_partitionId(partitionId),
      m_numBuckets(0),
      m_bucketViewLastBucketIdx(INVALID_INDEX),
      m_numEntities(0)
  {
    set_part_ordinals_from_host(partOrdinals);
  }

  KOKKOS_FUNCTION
  unsigned partition_id() const { return m_partitionId; }

  void set_partition_id(unsigned id) { m_partitionId = id; }

  KOKKOS_FUNCTION
  EntityRank get_rank() const { return m_rank; }

  // FIXME need to be adjusted when entities can be created on device (i.e. no entity info in bulkdata)
  bool add_entity(Entity entity)
  {
    auto bucketId = m_mesh->get_bulk_on_host().bucket(entity).bucket_id();
    auto bucketOrd = m_mesh->get_bulk_on_host().bucket_ordinal(entity);
    auto rank = m_mesh->get_bulk_on_host().entity_rank(entity);

    if (bucketId < m_deviceBucketRepo->m_numBuckets[rank]) {
      auto newBucket = m_deviceBucketRepo->construct_new_bucket(rank, superset_part_ordinals());
      add_bucket(newBucket);

      // TODO: Add entity to the bucket (newBucket) and update bucket connectivities
      return true;
    } else if (bucketOrd < m_deviceBucketRepo->m_buckets[rank](bucketId).size()) {

      // TODO: Add entity to the found bucket and update bucket connectivities
      return true;
    }

    return false;
  }

  void add_bucket(DeviceBucket* bucket)
  {
    bucket->m_owningPartitionId = m_partitionId;

    init_or_expand_bucket_view();

    DeviceBucketWrapper<NgpMemSpace> newBucketWrapper(bucket);

    auto newBucketIdx = get_next_avail_bucket_idx();
    m_buckets(newBucketIdx) = newBucketWrapper;

    update_partition_meta_bucket_added(bucket, newBucketIdx);
  }

  void remove_bucket(DeviceBucket* bucket)
  {
    m_deviceBucketRepo->invalidate_bucket(bucket);
    update_partition_meta_bucket_removed(bucket);
  }

  void sort_buckets()
  {
    using ExecSpace = typename NgpMemSpace::execution_space;
    DeviceBucketPtrView compactBucketPtrUView(m_buckets.data(), get_active_bucket_span());
    if (!Kokkos::Experimental::is_sorted(ExecSpace{}, compactBucketPtrUView)) {
      Kokkos::sort(ExecSpace{}, compactBucketPtrUView);
    }
  }

  void set_part_ordinals_from_host(std::vector<PartOrdinal> const& partOrdinals)
  {
    m_partOrdinals = PartOrdinalViewType<NgpMemSpace>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PartOrdinals"), partOrdinals.size());
    auto hostPartOrdinals = HostPartOrdinalViewType(partOrdinals.data(), partOrdinals.size());
    Kokkos::deep_copy(m_partOrdinals, hostPartOrdinals);
  }

  void init_bucket_view(unsigned numBuckets = initialDeviceBucketViewCapacity)
  {
    STK_ThrowRequireMsg(m_buckets.extent(0) == 0, "Bucket view in Device Partition is already initialized. Use init_or_expand_bucket_view() instead");
    m_buckets = DeviceBucketPtrView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketViewInPartition"), numBuckets);
  }

  void init_or_expand_bucket_view(unsigned numBuckets = initialDeviceBucketViewCapacity)
  {
    if (m_buckets.extent(0) == 0) {
      m_buckets = DeviceBucketPtrView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketViewInPartition"), numBuckets);
    } else if (m_numBuckets >= m_buckets.extent(0)) {
      Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_buckets, m_buckets.extent(0)*2);
    }
  }

  void reset_bucket_count()
  {
    if (has_no_buckets()) {
      m_partitionId = INVALID_PARTITION_ID;
    }
    m_numBuckets = 0;
    m_bucketViewLastBucketIdx = INVALID_INDEX;
  }

  void set_to_match_bucket_count()
  {
    if (has_no_buckets()) {
      reset_bucket_count();
    } else {
      Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_buckets, m_numBuckets);
    }
  }

  void reset_partition_id_in_owned_buckets()
  {
    for (unsigned i = 0; i < m_numBuckets; ++i) {
      auto bucket = m_buckets(i).bucketPtr;
      bucket->m_owningPartitionId = m_partitionId;
    }
  }

  void update_partition_meta_bucket_added(DeviceBucket* addedBucket, unsigned addedBucketIdx)
  {
    m_numEntities += addedBucket->size();
    m_numBuckets++;
    m_bucketViewLastBucketIdx = addedBucketIdx;
  }

  void update_partition_meta_bucket_removed(DeviceBucket* removedBucket)
  {
    m_numBuckets--;
    m_numEntities -= removedBucket->size();
    removedBucket->m_owningPartitionId = INVALID_PARTITION_ID;
  }

  KOKKOS_FUNCTION
  PartOrdinalViewType<NgpMemSpace> superset_part_ordinals() const { return m_partOrdinals; }

  KOKKOS_FUNCTION
  size_t num_buckets() const { return m_numBuckets; }

  KOKKOS_FUNCTION
  size_t num_entities() const { return m_numEntities; }

  KOKKOS_FUNCTION
  bool has_no_buckets() const { return m_numBuckets == 0; }

  KOKKOS_FUNCTION
  bool operator<(DevicePartition<NgpMemSpace> const& other) const {
    return m_partitionId < other.m_partitionId;
  }

  unsigned get_last_active_bucket_idx() const {
    return m_bucketViewLastBucketIdx;
  }

  unsigned get_active_bucket_span() const {
    auto lastBucketIdx = m_bucketViewLastBucketIdx;
    auto nextBucketIdx = (lastBucketIdx == INVALID_INDEX) ? 0 : lastBucketIdx+1;
    return nextBucketIdx;
  }

  unsigned get_next_avail_bucket_idx() const {
    auto lastBucketIdx = m_bucketViewLastBucketIdx;
    auto nextBucketIdx = (lastBucketIdx == INVALID_INDEX) ? 0 : lastBucketIdx+1;
    return nextBucketIdx;
  }

  void deep_copy_form_device_partition(DevicePartition<NgpMemSpace> const& rhs)
  {
    m_mesh = rhs.m_mesh;
    m_deviceBucketRepo = rhs.m_deviceBucketRepo;
    m_rank = rhs.m_rank;
    m_partitionId = rhs.m_partitionId;
    m_numBuckets = rhs.m_numBuckets;
    m_bucketViewLastBucketIdx = rhs.m_bucketViewLastBucketIdx;
    m_numEntities = rhs.m_numEntities;

    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_partOrdinals, rhs.m_partOrdinals.extent(0));
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_buckets, rhs.m_buckets.extent(0));

    auto execSpace = typename NgpMemSpace::execution_space{};
    Kokkos::Experimental::copy(execSpace, rhs.m_partOrdinals, m_partOrdinals);
    Kokkos::Experimental::copy(execSpace, rhs.m_buckets, m_buckets);
  }
  
  DeviceMeshT<NgpMemSpace>* m_mesh;
  DeviceBucketRepository<NgpMemSpace>* m_deviceBucketRepo;
  PartOrdinalViewType<NgpMemSpace> m_partOrdinals;
  DeviceBucketPtrView m_buckets;
  EntityRank m_rank;
  unsigned m_partitionId;
  unsigned m_numBuckets;
  unsigned m_bucketViewLastBucketIdx;
  unsigned m_numEntities;
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
    if (m_searchView(i).partition_id() == INVALID_PARTITION_ID) { return; }

    auto partOrdinals = m_searchView(i).superset_part_ordinals();

    if (m_partOrdinalView.extent(0) != partOrdinals.extent(0)) { return; }
  
    for (unsigned j = 0; j < m_partOrdinalView.extent(0); ++j) {
      if (partOrdinals(j) != m_partOrdinalView(j)) {
        return;
      }
    }
    update = m_searchView(i).partition_id();
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
