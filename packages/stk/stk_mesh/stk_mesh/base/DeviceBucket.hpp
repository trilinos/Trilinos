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

#ifndef STK_MESH_DEVICE_BUCKET_HPP
#define STK_MESH_DEVICE_BUCKET_HPP

#include <stk_util/stk_config.h>
#include <stk_util/util/StridedArray.hpp>
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/baseImpl/NgpMeshImpl.hpp"
#include "Kokkos_Core.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace mesh {

template <typename NgpMemSpace> class DeviceMeshT;

namespace impl {
template <typename NgpMemSpace>
class DevicePartition;

constexpr unsigned INVALID_INDEX = std::numeric_limits<unsigned>::max();

struct EntityCompareInvalidAtEnd
{
  using entity_value_type = stk::mesh::Entity::entity_value_type;
  static_assert(entity_value_type(Entity::InvalidEntity)-1 == std::numeric_limits<entity_value_type>::max(), 
                "InvalidEntity must wrap around to large value");

  KOKKOS_INLINE_FUNCTION
  constexpr bool operator()(const Entity& lhs, const Entity& rhs) const
  {
    return lhs.local_offset() - 1 < rhs.local_offset() - 1;
  }
};

}

template <typename BucketNgpMemSpace>
struct DeviceBucketT {
  using ConnectedNodes    = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedOrdinals = util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
  using Permutations      = util::StridedArray<const stk::mesh::Permutation>;

  KOKKOS_FUNCTION
  DeviceBucketT()
    : m_owningMesh(nullptr),
      m_owningPartitionId(INVALID_PARTITION_ID),
      m_bucketId(INVALID_BUCKET_ID),
      m_bucketSize(0),
      m_bucketCapacity(0),
      m_activeEntitySpan(0),
      m_bucketTopology(),
      m_entityRank(stk::topology::INVALID_RANK)
  {}

  KOKKOS_INLINE_FUNCTION
  unsigned bucket_id() const { return m_bucketId; }

  void set_bucket_id(unsigned id) { m_bucketId = id; }

  KOKKOS_INLINE_FUNCTION
  size_t size() const { return m_bucketSize; }

  KOKKOS_INLINE_FUNCTION
  size_t capacity() const { return m_bucketCapacity; }

  KOKKOS_INLINE_FUNCTION
  stk::mesh::EntityRank entity_rank() const { return m_entityRank; }

  KOKKOS_FUNCTION
  stk::topology topology() const { return m_bucketTopology; }

  KOKKOS_INLINE_FUNCTION
  unsigned partition_id() const { return m_owningPartitionId; }

  KOKKOS_INLINE_FUNCTION
  bool is_active() const { return m_bucketId != INVALID_BUCKET_ID; }

  KOKKOS_INLINE_FUNCTION
  unsigned get_active_entity_span() const { return m_activeEntitySpan; }

  KOKKOS_INLINE_FUNCTION
  ConnectedEntities get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  KOKKOS_INLINE_FUNCTION
  ConnectedOrdinals get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  KOKKOS_INLINE_FUNCTION
  Permutations get_connected_permutations(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const;

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_edges(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_faces(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_elements(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  stk::mesh::Entity operator[](unsigned offsetIntoBucket) const {
    STK_NGP_ThrowAssert(offsetIntoBucket < m_entities.size());
    return m_entities(offsetIntoBucket);
  }

  KOKKOS_FUNCTION
  bool member(stk::mesh::PartOrdinal partOrdinal) const
  {
    for(unsigned i=0; i<m_partOrdinals.size(); i++) {
      if(m_partOrdinals(i) == partOrdinal) {
        return true;
      }
    }
    return false;
  }

  KOKKOS_FUNCTION
  const Kokkos::pair<const stk::mesh::PartOrdinal*,const stk::mesh::PartOrdinal*> superset_part_ordinals() const
  {
    STK_NGP_ThrowAssert(m_partOrdinals.size() > 0);
    return Kokkos::pair<const stk::mesh::PartOrdinal*,const stk::mesh::PartOrdinal*>(
             m_partOrdinals.data(), m_partOrdinals.data()+m_partOrdinals.size()
           );
  }

  KOKKOS_FUNCTION
  bool operator<(DeviceBucketT<BucketNgpMemSpace> const& rhs) const
  {
    if (m_bucketId == INVALID_BUCKET_ID) {
      return false;
    } else if (rhs.m_bucketId == INVALID_BUCKET_ID) {
      return true;
    } else {
      return impl::DevicePartOrdinalLess{}(get_part_ordinals(), rhs.get_part_ordinals());
    }
  }

  void initialize_bucket_attributes(const stk::mesh::Bucket &bucket);
  void initialize_part_ordinals_from_host(const stk::mesh::Bucket &bucket);
  void update_entity_data_from_host(const stk::mesh::Bucket &bucket);
  void update_sparse_connectivity_from_host(const stk::mesh::Bucket &bucket);

  void resize_device_views(const stk::mesh::Bucket &bucket);
  std::pair<unsigned, unsigned> scan_entities_for_nodal_connectivity(const stk::mesh::Bucket & bucket);

  void clone_device_bucket(DeviceBucketT<BucketNgpMemSpace> const& rhs);

  void sort_entities();

  void init_entity_view();

  void add_entity(Entity entity);

  KOKKOS_FUNCTION
  void remove_entity(Entity entity);

  KOKKOS_FUNCTION
  void remove_entity(unsigned bucketOrd);

  // FIXME. Entity view could be sparse during a mesh mod.
  bool is_full() const { return get_next_avail_entity_idx() >= capacity(); }

  KOKKOS_FUNCTION
  const PartOrdinalViewType<BucketNgpMemSpace>& get_part_ordinals() const { return m_partOrdinals; }

  void set_part_ordinals(PartOrdinalViewType<BucketNgpMemSpace>& view) { m_partOrdinals = view; }

  void update_bucket_meta_entity_added()
  {
    m_bucketSize++;
    m_activeEntitySpan++;
  }

  KOKKOS_FUNCTION
  void update_bucket_meta_entity_removed()
  {
    Kokkos::atomic_dec(&m_bucketSize);
  }

  KOKKOS_INLINE_FUNCTION
  void update_bucket_meta_set_entity_span_to_active_count()
  {
    m_activeEntitySpan = m_bucketSize;
  }

  unsigned get_next_avail_entity_idx() const {
    auto bucketSpan = m_activeEntitySpan;
    auto nextEntityIdx = (bucketSpan >= capacity()) ? impl::INVALID_INDEX : bucketSpan;
    return nextEntityIdx;
  }

  EntityViewType<BucketNgpMemSpace> m_entities;
  BucketConnectivityType<BucketNgpMemSpace> m_nodeConnectivity;
  OrdinalViewType<BucketNgpMemSpace> m_nodeConnectivityOffsets;
  OrdinalViewType<BucketNgpMemSpace> m_nodeOrdinals;
  Unsigned2dViewType<BucketNgpMemSpace> m_sparseConnectivityOffsets;
  BucketConnectivityType<BucketNgpMemSpace> m_sparseConnectivity;
  OrdinalViewType<BucketNgpMemSpace> m_sparseConnectivityOrdinals;
  PermutationViewType<BucketNgpMemSpace> m_sparseConnectivityPermutations;
  PartOrdinalViewType<BucketNgpMemSpace> m_partOrdinals;
  
  const stk::mesh::DeviceMeshT<BucketNgpMemSpace>* m_owningMesh;
  unsigned m_owningPartitionId;
  unsigned m_bucketId;
  unsigned m_bucketSize;
  unsigned m_bucketCapacity;
  unsigned m_activeEntitySpan;
  stk::topology m_bucketTopology;
  stk::mesh::EntityRank m_entityRank;
};

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::add_entity(Entity entity)
{
  STK_ThrowAssertMsg(!is_full(), "DeviceBucket must not be full to add an entity");
  auto nextIdx = get_next_avail_entity_idx();
  STK_ThrowAssertMsg(nextIdx != impl::INVALID_INDEX, "Bucket is full");

  Kokkos::parallel_for(1, KOKKOS_CLASS_LAMBDA(const int) {
    m_entities(nextIdx) = entity;
  });
  Kokkos::fence();

  update_bucket_meta_entity_added();

  // Update connectivities based on added entity?
}

template <typename BucketNgpMemSpace>
KOKKOS_FUNCTION
void DeviceBucketT<BucketNgpMemSpace>::remove_entity(Entity entity)
{
  for (unsigned i = 0; i < m_activeEntitySpan; ++i) {
    if (m_entities(i) == entity) {
      m_entities(i) = Entity();
    }
  }
  update_bucket_meta_entity_removed();
}

template <typename BucketNgpMemSpace>
KOKKOS_FUNCTION
void DeviceBucketT<BucketNgpMemSpace>::remove_entity(unsigned bucketOrd)
{
  m_entities(bucketOrd) = Entity();
  update_bucket_meta_entity_removed();
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::sort_entities()
{
  using EntityUViewType = Kokkos::View<Entity*, BucketNgpMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  EntityUViewType compactEntityUView(m_entities.data(), get_active_entity_span());
  if (!Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactEntityUView, impl::EntityCompareInvalidAtEnd{})) {
    Kokkos::sort(compactEntityUView, impl::EntityCompareInvalidAtEnd{});
  }
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactEntityUView, impl::EntityCompareInvalidAtEnd{}));

  update_bucket_meta_set_entity_span_to_active_count();
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::init_entity_view()
{
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_entities, capacity());
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::clone_device_bucket(DeviceBucketT<BucketNgpMemSpace> const& rhs)
{
  m_owningMesh = rhs.m_owningMesh;
  m_owningPartitionId = rhs.m_owningPartitionId;
  m_bucketId = rhs.m_bucketId;
  m_bucketSize = rhs.m_bucketSize;
  m_bucketCapacity = rhs.m_bucketCapacity;
  m_bucketTopology = rhs.m_bucketTopology;
  m_entityRank = rhs.m_entityRank;

  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_entities, rhs.m_entities.extent(0));
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_nodeConnectivity, rhs.m_nodeConnectivity.extent(0));
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_nodeConnectivityOffsets, rhs.m_nodeConnectivityOffsets.extent(0));
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_nodeOrdinals, rhs.m_nodeOrdinals.extent(0));
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_sparseConnectivityOffsets, rhs.m_sparseConnectivityOffsets.extent(0), rhs.m_sparseConnectivityOffsets.extent(1));
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_sparseConnectivity, rhs.m_sparseConnectivity.extent(0));
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_sparseConnectivityOrdinals, rhs.m_sparseConnectivityOrdinals.extent(0));
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_partOrdinals, rhs.m_partOrdinals.extent(0));

  auto execSpace = typename BucketNgpMemSpace::execution_space{};
  Kokkos::Experimental::copy(execSpace, rhs.m_entities, m_entities);
  Kokkos::Experimental::copy(execSpace, rhs.m_nodeConnectivity, m_nodeConnectivity);
  Kokkos::Experimental::copy(execSpace, rhs.m_nodeConnectivityOffsets, m_nodeConnectivityOffsets);
  Kokkos::Experimental::copy(execSpace, rhs.m_nodeOrdinals, m_nodeOrdinals);
  Kokkos::deep_copy(m_sparseConnectivityOffsets, rhs.m_sparseConnectivityOffsets);
  Kokkos::Experimental::copy(execSpace, rhs.m_sparseConnectivity, m_sparseConnectivity);
  Kokkos::Experimental::copy(execSpace, rhs.m_sparseConnectivityOrdinals, m_sparseConnectivityOrdinals);
  Kokkos::Experimental::copy(execSpace, rhs.m_partOrdinals, m_partOrdinals);
}

template <typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::ConnectedEntities
DeviceBucketT<BucketNgpMemSpace>::get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  if (connectedRank == stk::topology::NODE_RANK) {
    const unsigned numNodes = m_nodeConnectivityOffsets(offsetIntoBucket+1)-m_nodeConnectivityOffsets(offsetIntoBucket);
    const size_t nodeOffset = m_nodeConnectivityOffsets(offsetIntoBucket);
    return ConnectedEntities(&m_nodeConnectivity(nodeOffset), numNodes, 1);
  }

  const unsigned offset = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket);
  const unsigned length = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket+1) - offset;
  return ConnectedEntities(&m_sparseConnectivity(offset), length, 1);
}

template <typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::ConnectedOrdinals
DeviceBucketT<BucketNgpMemSpace>::get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  if (connectedRank == stk::topology::NODE_RANK) {
    const unsigned numNodes = m_nodeConnectivityOffsets(offsetIntoBucket+1)-m_nodeConnectivityOffsets(offsetIntoBucket);
    return ConnectedOrdinals(m_nodeOrdinals.data(), numNodes, 1);
  }

  const unsigned offset = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket);
  const unsigned length = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket+1) - offset;
  return ConnectedOrdinals(&m_sparseConnectivityOrdinals(offset), length, 1);
}

template <typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::Permutations
DeviceBucketT<BucketNgpMemSpace>::get_connected_permutations(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  const unsigned offset = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket);
  const unsigned length = m_sparseConnectivityOffsets(connectedRank,offsetIntoBucket+1) - offset;
  if (m_sparseConnectivityPermutations.size() <= offset) {
    return Permutations(nullptr, 0);
  }

  return Permutations(&m_sparseConnectivityPermutations(offset), length, 1);
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::update_sparse_connectivity_from_host(const stk::mesh::Bucket &bucket)
{
  Kokkos::Profiling::pushRegion("update_sparse_connectivity_from_host()");

  typename Unsigned2dViewType<BucketNgpMemSpace>::host_mirror_type hostConnectivityOffsets("hostConnectivityOffsets", 0,0); 
  Kokkos::resize(Kokkos::WithoutInitializing, hostConnectivityOffsets, stk::topology::NUM_RANKS, bucket.size()+1);
  Kokkos::resize(Kokkos::WithoutInitializing, m_sparseConnectivityOffsets, stk::topology::NUM_RANKS, bucket.size()+1);
  typename BucketConnectivityType<BucketNgpMemSpace>::host_mirror_type hostConnectivity("hostConnectivity", 0);
  typename OrdinalViewType<BucketNgpMemSpace>::host_mirror_type hostConnectivityOrdinals("hostConnectivityOrdinals", 0);
  typename PermutationViewType<BucketNgpMemSpace>::host_mirror_type hostConnectivityPermutations("hostConnectivityPermutations", 0);

  const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bucket.mesh().mesh_meta_data().entity_rank_count());

  unsigned offset = 0;
  for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; ++connectedRank) {
    for(unsigned i=0; i<bucket.size(); ++i) {
      hostConnectivityOffsets(connectedRank,i) = offset;
      offset += bucket.num_connectivity(i, connectedRank);
    }
    hostConnectivityOffsets(connectedRank,bucket.size()) = offset;
  }

  Kokkos::resize(Kokkos::WithoutInitializing, hostConnectivity, offset);
  Kokkos::resize(Kokkos::WithoutInitializing, m_sparseConnectivity, offset);
  Kokkos::resize(Kokkos::WithoutInitializing, hostConnectivityOrdinals, offset);
  Kokkos::resize(Kokkos::WithoutInitializing, m_sparseConnectivityOrdinals, offset);
  const bool hasPermutations = bucket.has_permutation(stk::topology::EDGE_RANK)
                            || bucket.has_permutation(stk::topology::FACE_RANK)
                            || bucket.has_permutation(stk::topology::ELEM_RANK);
  if (hasPermutations) {
    Kokkos::resize(Kokkos::WithoutInitializing, hostConnectivityPermutations, offset);
    Kokkos::resize(Kokkos::WithoutInitializing, m_sparseConnectivityPermutations, offset);
  }

  offset = 0;
  for(stk::mesh::EntityRank connectedRank=stk::topology::EDGE_RANK; connectedRank<endRank; ++connectedRank) {
    for(unsigned i=0; i<bucket.size(); ++i) {
      const unsigned numConn = bucket.num_connectivity(i, connectedRank);
      const Entity* conn = bucket.begin(i, connectedRank);
      const ConnectivityOrdinal* ords = bucket.begin_ordinals(i, connectedRank);
      const Permutation* perms = bucket.begin_permutations(i, connectedRank);

      for(unsigned nc=0; nc<numConn; ++nc) {
        hostConnectivity(offset+nc) = conn[nc];
        hostConnectivityOrdinals(offset+nc) = ords[nc];
        if (hasPermutations) {
          hostConnectivityPermutations(offset+nc) = perms[nc];
        }
      }
      offset += numConn;
    }
  }

  Kokkos::deep_copy(m_sparseConnectivityOffsets, hostConnectivityOffsets);
  Kokkos::deep_copy(m_sparseConnectivity, hostConnectivity);
  Kokkos::deep_copy(m_sparseConnectivityOrdinals, hostConnectivityOrdinals);
  if (hasPermutations) {
    Kokkos::deep_copy(m_sparseConnectivityPermutations, hostConnectivityPermutations);
  }

  Kokkos::Profiling::popRegion();
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::initialize_bucket_attributes(const stk::mesh::Bucket &bucket)
{
  m_bucketId = bucket.bucket_id();
  m_bucketCapacity = bucket.capacity();
  m_bucketSize = bucket.size();
  m_entityRank = bucket.entity_rank();
  m_bucketTopology = bucket.topology();
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::initialize_part_ordinals_from_host(const stk::mesh::Bucket &bucket)
{
  const stk::mesh::PartVector& parts = bucket.supersets();
  m_partOrdinals = PartOrdinalViewType<BucketNgpMemSpace>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "PartOrdinals"),
                                                          parts.size());
  auto hostPartOrdinals = HostPartOrdinalViewType(bucket.superset_part_ordinals().first, parts.size());
  Kokkos::deep_copy(m_partOrdinals, hostPartOrdinals);
}

template <typename BucketNgpMemSpace>
std::pair<unsigned, unsigned>
DeviceBucketT<BucketNgpMemSpace>::scan_entities_for_nodal_connectivity(const stk::mesh::Bucket & bucket)
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

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::resize_device_views(const stk::mesh::Bucket & bucket)
{
  Kokkos::Profiling::pushRegion("resize_device_views()");

  Kokkos::Profiling::pushRegion("set node ordinals");

  const auto [maxNodesPerEntity, totalNumConnectedNodes] = scan_entities_for_nodal_connectivity(bucket);

  if (m_nodeOrdinals.size() != maxNodesPerEntity) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_nodeOrdinals, static_cast<size_t>(maxNodesPerEntity));
    OrdinalViewType<BucketNgpMemSpace>& nodeOrds = m_nodeOrdinals; //local var to avoid implicit this capture
    Kokkos::parallel_for(Kokkos::RangePolicy<stk::ngp::ExecSpace>(0, maxNodesPerEntity),
      KOKKOS_LAMBDA(const int i) {
        nodeOrds(i) = static_cast<stk::mesh::ConnectivityOrdinal>(i);
      });
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("bucket entities");
  if (m_entities.size() != m_bucketCapacity) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_entities, m_bucketCapacity);
    STK_ThrowRequireMsg(m_bucketCapacity > 0, "bucket capacity must be greater than 0");
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("nodal connectivity");
  if (m_nodeConnectivity.size() != totalNumConnectedNodes) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_nodeConnectivity, totalNumConnectedNodes);
  }

  if (m_nodeConnectivityOffsets.size() != m_bucketCapacity+1) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_nodeConnectivityOffsets, m_bucketCapacity+1);
  }
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::update_entity_data_from_host(const stk::mesh::Bucket &bucket)
{
  Kokkos::Profiling::pushRegion("update_entity_data_from_host()");

  m_bucketSize = bucket.size();
  m_bucketCapacity = bucket.capacity();
  m_activeEntitySpan = bucket.size();

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

} }

#endif
