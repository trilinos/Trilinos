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

namespace stk {
namespace mesh {

template<typename NgpMemSpace>
class DeviceMeshT;

namespace impl {

constexpr size_t initialBucketViewCapacity = 32;
constexpr size_t initialPartitionViewCapacity = 32;

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

  KOKKOS_DEFAULTED_FUNCTION
  DevicePartition() = default;

  DevicePartition(DeviceMeshT<NgpMemSpace>* deviceMesh, DeviceBucketRepository<NgpMemSpace>* deviceBucketRepo,
                  PartOrdinalViewType<NgpMemSpace> const& partOrdinals, EntityRank rank, unsigned partitionId)
    : m_mesh(deviceMesh),
      m_deviceBucketRepo(deviceBucketRepo),
      m_partOrdinals(partOrdinals),
      m_rank(rank),
      m_partitionId(partitionId),
      m_numBuckets(0),
      m_bucketViewEndIdx(0),
      m_numEntities(0)
  {
    init_or_resize_bucket_view();
  }
  
  KOKKOS_DEFAULTED_FUNCTION DevicePartition(DevicePartition const&) = default;
  KOKKOS_DEFAULTED_FUNCTION DevicePartition(DevicePartition&&) = default;
  KOKKOS_DEFAULTED_FUNCTION DevicePartition& operator=(DevicePartition&) = default;
  KOKKOS_DEFAULTED_FUNCTION DevicePartition& operator=(DevicePartition const&) = default;

  KOKKOS_FUNCTION
  unsigned partition_id() const { return m_partitionId; }

  // FIXME need to be adjusted when entities can be created on device (i.e. no entity info in bulkdata)
  bool add_entity(Entity entity)
  {
    auto bucketId = m_mesh->get_bulk_on_host().bucket(entity).bucket_id();
    auto bucketOrd = m_mesh->get_bulk_on_host().bucket_ordinal(entity);
    auto rank = m_mesh->get_bulk_on_host().entity_rank(entity);

    if (bucketId < m_deviceBucketRepo->m_numBuckets[rank]) {
      auto newBucket = m_deviceBucketRepo->allocate_bucket(rank, superset_part_ordinals());
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

    init_or_resize_bucket_view();

    DeviceBucketWrapper<NgpMemSpace> newBucketWrapper(bucket);
    m_buckets(m_numBuckets) = newBucketWrapper;

    m_numEntities += bucket->size();
    m_numBuckets++;
    m_bucketViewEndIdx++;
  }

  void remove_bucket(DeviceBucket* bucket)
  {
    m_deviceBucketRepo->deallocate_bucket(bucket);
    m_numBuckets--;
    m_numEntities -= bucket->size();
    bucket->m_owningPartitionId = INVALID_PARTITION_ID;
  }

  void sort_buckets()
  {
    using ExecSpace = typename NgpMemSpace::execution_space;
    DeviceBucketPtrView compactBucketPtrUView(m_buckets.data(), m_bucketViewEndIdx - 1);
    if (!Kokkos::Experimental::is_sorted(ExecSpace{}, compactBucketPtrUView)) {
      Kokkos::sort(ExecSpace{}, compactBucketPtrUView);
    }
  }

  void init_or_resize_bucket_view()
  {
    if (m_buckets.extent(0) == 0) {
      m_buckets = DeviceBucketPtrView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceBucketViewInPartition"), initialBucketViewCapacity);
    } else if (m_numBuckets >= m_buckets.extent(0)) {
      Kokkos::resize(m_buckets, m_buckets.extent(0)*2);
    }
  }

  void reset_bucket_count()
  {
    if (has_no_buckets()) {
      m_partitionId = INVALID_PARTITION_ID;
    }
    m_numBuckets = 0;
    m_bucketViewEndIdx = 0;
  }

  void compact_bucket_view_and_reset_bucket_count()
  {
    if (has_no_buckets()) {
      m_partitionId = INVALID_PARTITION_ID;
    } else {
      Kokkos::resize(Kokkos::WithoutInitializing, m_buckets, m_numBuckets);
      m_numBuckets = 0;
      m_bucketViewEndIdx = 0;
    }
  }

  void reset_partition_id_in_owned_buckets()
  {
    for (unsigned i = 0; i < m_numBuckets; ++i) {
      auto bucket = m_buckets(i).bucketPtr;
      bucket->m_owningPartitionId = m_partitionId;
    }
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
  
  DeviceMeshT<NgpMemSpace>* m_mesh;
  DeviceBucketRepository<NgpMemSpace>* m_deviceBucketRepo;
  PartOrdinalViewType<NgpMemSpace> m_partOrdinals;
  DeviceBucketPtrView m_buckets;
  EntityRank m_rank;
  unsigned m_partitionId;
  size_t m_numBuckets;
  size_t m_bucketViewEndIdx;
  size_t m_numEntities;
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
