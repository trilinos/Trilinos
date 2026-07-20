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

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/base/BucketConnectivity.hpp"
#include "stk_mesh/baseImpl/NgpMeshImpl.hpp"
#include "stk_mesh/baseImpl/MeshConnectivity.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/stk_config.h"
#include "stk_util/util/StridedArray.hpp"
#include "Kokkos_Core.hpp"

namespace stk {
namespace mesh {

template <typename NgpMemSpace> class DeviceMeshT;

namespace impl {
template <typename NgpMemSpace>
class DevicePartition;

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

struct EntityAndIndex
{
  KOKKOS_INLINE_FUNCTION
  EntityAndIndex() :
    entity{},
    idx{static_cast<unsigned>(-1)}
  {}

  KOKKOS_INLINE_FUNCTION
  EntityAndIndex(Entity entityArg, unsigned idxArg) :
    entity(entityArg),
    idx(idxArg)
  {}

  KOKKOS_INLINE_FUNCTION
  operator Entity() const { return entity; }

  Entity entity;
  unsigned idx;
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
      m_meshConn(),
      m_owningPartitionId(INVALID_PARTITION_ID),
      m_bucketId(INVALID_BUCKET_ID),
      m_bucketSize(0),
      m_bucketCapacity(0),
      m_activeEntitySpan(0),
      m_bucketTopology(),
      m_entityRank(stk::topology::INVALID_RANK),
      m_entityEndRank(stk::topology::INVALID_RANK),
      m_isModified(false),
      m_hasPermutations(false)
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
  void set_inactive() { m_bucketId = INVALID_BUCKET_ID; }

  KOKKOS_INLINE_FUNCTION
  bool is_modified() const { return m_isModified; }

  KOKKOS_INLINE_FUNCTION
  void set_modified() { m_isModified = true; }

  KOKKOS_INLINE_FUNCTION
  void clear_modified() { m_isModified = false; }

  KOKKOS_INLINE_FUNCTION
  void set_has_permutations(bool hasPermutations) { m_hasPermutations = hasPermutations; }

  KOKKOS_INLINE_FUNCTION
  bool has_permutations() const { return m_hasPermutations; }

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

  void sync_to_host(stk::mesh::Bucket& bucket);

  void clone_device_bucket(DeviceBucketT<BucketNgpMemSpace> const& rhs);

  void sort_entities();

  void init_entity_view();

  void add_entity(Entity entity);

  KOKKOS_FUNCTION
  void clear_entities();

  KOKKOS_FUNCTION
  void add_entity(Entity entity, unsigned bucketOrd);

  KOKKOS_FUNCTION
  void remove_entity(Entity entity);

  KOKKOS_FUNCTION
  void remove_entity(unsigned bucketOrd);

  // FIXME. Entity view could be sparse during a mesh mod.
  KOKKOS_INLINE_FUNCTION
  bool is_full() const { return get_next_avail_entity_idx() >= capacity(); }

  KOKKOS_FUNCTION
  const PartOrdinalViewType<BucketNgpMemSpace>& get_part_ordinals() const { return m_partOrdinals; }

  void set_part_ordinals(PartOrdinalViewType<BucketNgpMemSpace>& view) { m_partOrdinals = view; }

  KOKKOS_FUNCTION
  void update_bucket_meta_entity_added()
  {
    Kokkos::atomic_inc(&m_bucketSize);
    Kokkos::atomic_inc(&m_activeEntitySpan);
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

  KOKKOS_INLINE_FUNCTION
  unsigned get_next_avail_entity_idx() const {
    auto bucketSpan = m_activeEntitySpan;
    auto nextEntityIdx = (bucketSpan >= capacity()) ? impl::INVALID_INDEX : bucketSpan;
    return nextEntityIdx;
  }

  KOKKOS_INLINE_FUNCTION
  void set_end_rank(EntityRank rank) { m_entityEndRank = rank; }

  KOKKOS_INLINE_FUNCTION
  EntityRank get_end_rank() const { return m_entityEndRank; }

  void set_mesh_connectivity(impl::MeshConnectivity<stk::ngp::UVMDeviceSpace>& meshConn)
  {
    m_meshConn = meshConn;
  }

  const DeviceMeshT<BucketNgpMemSpace>* m_owningMesh;
  impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> m_meshConn;
  EntityViewType<BucketNgpMemSpace> m_entities;
  PartOrdinalViewType<BucketNgpMemSpace> m_partOrdinals;

  unsigned m_owningPartitionId;
  unsigned m_bucketId;
  unsigned m_bucketSize;
  unsigned m_bucketCapacity;
  unsigned m_activeEntitySpan;
  stk::topology m_bucketTopology;
  EntityRank m_entityRank;
  EntityRank m_entityEndRank;
  bool m_isModified;
  bool m_hasPermutations;

  void permute_field_data(const Kokkos::View<const impl::EntityAndIndex*, BucketNgpMemSpace>& sortedEntities);

  private:
    void update_entity_data_from_device(stk::mesh::Bucket& bucket);
    void update_sparse_connectivity_from_device(stk::mesh::Bucket &bucket);
};

template <typename BucketNgpMemSpace>
KOKKOS_FUNCTION
void DeviceBucketT<BucketNgpMemSpace>::clear_entities()
{
  m_bucketSize = 0;
  update_bucket_meta_set_entity_span_to_active_count();
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::add_entity(Entity entity)
{
  auto nextIdx = get_next_avail_entity_idx();
  KOKKOS_IF_ON_HOST((
    STK_ThrowAssertMsg(!is_full(), "DeviceBucket must not be full to add an entity");
    STK_ThrowAssertMsg(nextIdx != impl::INVALID_INDEX, "Bucket is full");
  ))

  add_entity(entity, nextIdx);
}

template <typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
void DeviceBucketT<BucketNgpMemSpace>::add_entity(Entity entity, unsigned destBucketOrd)
{
  m_entities(destBucketOrd) = entity;

  update_bucket_meta_entity_added();

  // Update connectivities based on added entity?
}

template <typename BucketNgpMemSpace>
KOKKOS_FUNCTION
void DeviceBucketT<BucketNgpMemSpace>::remove_entity(Entity entity)
{
  for (unsigned i = 0; i < m_activeEntitySpan; ++i) {
    if (m_entities(i) == entity) {
      m_entities(i) = Entity{};
    }
  }
  update_bucket_meta_entity_removed();
}

template <typename BucketNgpMemSpace>
KOKKOS_FUNCTION
void DeviceBucketT<BucketNgpMemSpace>::remove_entity(unsigned bucketOrd)
{
  m_entities(bucketOrd) = Entity{};
  update_bucket_meta_entity_removed();
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::sort_entities()
{
  using EntityUViewType = Kokkos::View<Entity*, BucketNgpMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  EntityUViewType compactEntityUView(m_entities.data(), get_active_entity_span());
  if (!Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactEntityUView, impl::EntityCompareInvalidAtEnd{})) {
    Kokkos::View<impl::EntityAndIndex*, BucketNgpMemSpace> entitiesToSort(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entities_to_sort"), compactEntityUView.size());

    Kokkos::parallel_for("copy_entities", entitiesToSort.size(),
      KOKKOS_CLASS_LAMBDA (const unsigned idx)
      {
        entitiesToSort(idx) = impl::EntityAndIndex{m_entities[idx], idx};
      }
    );

    Kokkos::sort(entitiesToSort, impl::EntityCompareInvalidAtEnd{});

    Kokkos::deep_copy(compactEntityUView, entitiesToSort);

    permute_field_data(entitiesToSort);
  }
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{}, compactEntityUView, impl::EntityCompareInvalidAtEnd{}));

  update_bucket_meta_set_entity_span_to_active_count();
}


template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::permute_field_data(const Kokkos::View<const impl::EntityAndIndex*, BucketNgpMemSpace>& sortedEntities)
{
  using ExecSpace = typename BucketNgpMemSpace::execution_space;
  using TeamPolicy = stk::ngp::TeamPolicy<ExecSpace>;
  using TeamMember = typename TeamPolicy::member_type;
  using FieldDataBytesType = FieldDataBytes<ngp::DeviceSpace>;
  using ScratchView = Kokkos::View<std::byte**, typename TeamMember::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  const BulkData& bulk = m_owningMesh->get_bulk_on_host();
  const FieldVector& fields = bulk.mesh_meta_data().get_fields(entity_rank());
  Kokkos::View<FieldDataBytesType*, BucketNgpMemSpace> fieldDataBytes(Kokkos::view_alloc(Kokkos::WithoutInitializing, "field_data_bytes"), fields.size());
  assemble_field_data_bytes_on_device<ngp::DeviceSpace>(fields, fieldDataBytes);

  unsigned maxBytesPerEntity = 0;
  for (unsigned i=0; i < fields.size(); ++i)
  {
    unsigned bytesPerEntity = field_scalars_per_entity(*(fields[i]), bucket_id()) * fields[i]->data_traits().size_of;
    maxBytesPerEntity = std::max(maxBytesPerEntity, bytesPerEntity);
  }

  TeamPolicy policy(fields.size(), Kokkos::AUTO);
  policy.set_scratch_size(0, Kokkos::PerTeam(maxBytesPerEntity*m_entities.size()));

  unsigned bucketId = bucket_id();
  auto func = KOKKOS_LAMBDA(const TeamMember& team)
  {
    auto fieldIdx = team.league_rank();

    FieldDataBytesType& fieldData = fieldDataBytes[fieldIdx];
    auto bucketBytes = fieldData.bucket_bytes<Layout::Left>(bucketId);
    if (bucketBytes.num_bytes() == 0)
    {
      return;
    }

    ScratchView valuesCopy(team.team_scratch(0), sortedEntities.size(), bucketBytes.num_bytes());
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, sortedEntities.size()),
      [&](const int entityIdx)
      {
        unsigned oldEntityIndex = sortedEntities(entityIdx).idx;
        for (auto byteIdx : bucketBytes.bytes())
        {
          valuesCopy(entityIdx, byteIdx()) = bucketBytes(static_cast<EntityIdx>(oldEntityIndex), byteIdx);
        }
      }
    );

    team.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, sortedEntities.size()),
      [&](const int entityIdx)
      {
        for (auto byteIdx : bucketBytes.bytes())
          bucketBytes(static_cast<EntityIdx>(entityIdx), byteIdx) = valuesCopy(entityIdx, byteIdx());
      }
    );
  };

  Kokkos::parallel_for("permute_bucket_field_data", policy, func);

  Kokkos::fence();
}


template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::init_entity_view()
{
  Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), m_entities, capacity());
}

template <typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::ConnectedEntities
DeviceBucketT<BucketNgpMemSpace>::get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  auto connEntsNew = m_meshConn.get_connected_entities(m_entities(offsetIntoBucket), connectedRank);
  return connEntsNew;
}

template <typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::ConnectedOrdinals
DeviceBucketT<BucketNgpMemSpace>::get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  auto connOrdsNew = m_meshConn.get_connected_ordinals(m_entities(offsetIntoBucket), connectedRank);
  return connOrdsNew;
}

template <typename BucketNgpMemSpace>
KOKKOS_INLINE_FUNCTION
typename DeviceBucketT<BucketNgpMemSpace>::Permutations
DeviceBucketT<BucketNgpMemSpace>::get_connected_permutations(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
  STK_NGP_ThrowAssert(connectedRank < stk::topology::NUM_RANKS);
  auto connPermsNew = m_meshConn.get_connected_permutations(m_entities(offsetIntoBucket), connectedRank);
  if (connPermsNew.size() == 0) {
    return Permutations(nullptr, 0);
  }

  return connPermsNew;
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::sync_to_host(stk::mesh::Bucket& bucket)
{
  update_entity_data_from_device(bucket);
  update_sparse_connectivity_from_device(bucket);
}

template <typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::initialize_bucket_attributes(const stk::mesh::Bucket &bucket)
{
  m_bucketId = bucket.bucket_id();
  m_bucketCapacity = bucket.capacity();
  m_bucketSize = bucket.size();
  m_entityRank = bucket.entity_rank();
  m_entityEndRank = static_cast<EntityRank>(bucket.mesh().mesh_meta_data().entity_rank_count());
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
void DeviceBucketT<BucketNgpMemSpace>::update_entity_data_from_host(const stk::mesh::Bucket &bucket)
{
  Kokkos::Profiling::pushRegion("update_entity_data_from_host()");

  m_bucketSize = bucket.size();
  m_bucketCapacity = bucket.capacity();
  m_activeEntitySpan = bucket.size();

  if (m_entities.size() != m_bucketCapacity) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_entities, m_bucketCapacity);
    STK_ThrowRequireMsg(m_bucketCapacity > 0, "bucket capacity must be greater than 0");
  }

  auto hostEntities = HostEntityViewType(bucket.begin(), m_bucketCapacity);
  Kokkos::deep_copy(m_entities, hostEntities);

  Kokkos::Profiling::popRegion();
}

template<typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::update_entity_data_from_device(stk::mesh::Bucket &bucket)
{
  Kokkos::Profiling::pushRegion("update_entity_data_from_device()");

  auto hostEntities = Kokkos::create_mirror_view(m_entities);
  Kokkos::deep_copy(hostEntities, m_entities);

  bucket.set_ngp_field_bucket_id(bucket_id());
  bucket.set_ngp_mesh_bucket_id(bucket_id());
  int size_difference = std::abs(static_cast<int>(bucket.size()) - static_cast<int>(size()));

  if (bucket.size() > size())
  {
    for (int i=0; i < size_difference; ++i)
    {
      bucket.remove_entity();
    }
  } else
  {
    while(bucket.capacity() < size()) {
      bucket.grow_capacity();
    }
    for (int i=0; i < size_difference; ++i)
    {
      bucket.add_entity();
    }
  }

  for (size_t i=0; i < size(); ++i)
  {
    bucket.m_entities[i] = hostEntities(i);
  }

  Kokkos::Profiling::popRegion();
}

template<typename BucketNgpMemSpace>
void DeviceBucketT<BucketNgpMemSpace>::update_sparse_connectivity_from_device(stk::mesh::Bucket &bucket)
{

  STK_ThrowAssertMsg(get_end_rank() == bucket.mesh().mesh_meta_data().entity_rank_count(),
                     "Device bucket end rank and host bucket end rank do not match");

  for (stk::mesh::EntityRank rank=stk::topology::BEGIN_RANK; rank < bucket.mesh().mesh_meta_data().entity_rank_count(); ++rank)
  {
    if (rank == entity_rank())
    {
      continue;
    }

    for (size_t i=0; i < size(); ++i)
    {
      auto connEntities = m_meshConn.get_connected_entities(bucket[i], rank);
      auto connOrdinals = m_meshConn.get_connected_ordinals(bucket[i], rank);
      auto connPerms = m_meshConn.get_connected_permutations(bucket[i], rank);

      bucket.replace_relations(i, rank, connEntities.size(), connEntities.data(), connOrdinals.data(), connPerms.data());
    }
  }
}

}
}

#endif
