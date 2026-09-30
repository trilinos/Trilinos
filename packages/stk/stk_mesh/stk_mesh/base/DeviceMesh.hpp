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

#ifndef STK_MESH_DEVICEMESH_HPP
#define STK_MESH_DEVICEMESH_HPP

#include <stk_util/stk_config.h>
#include "Kokkos_Macros.hpp"
#include "stk_mesh/base/DeviceFieldDataManagerBase.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include "Kokkos_Core.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/DeviceFieldDataManager.hpp"
#include <string>

#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/baseImpl/NgpMeshImpl.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/baseImpl/MeshConnectivity.hpp"
#include "stk_mesh/baseImpl/MeshConnUtils.hpp"
#include "stk_mesh/baseImpl/PartVectorUtils.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/util/StkNgpVector.hpp"
#include "stk_util/util/ReportHandler.hpp"

#include "stk_mesh/baseImpl/DeviceMeshViewVector.hpp"
#include "stk_mesh/baseImpl/Partition.hpp"
#include "stk_mesh/baseImpl/NgpMeshHostData.hpp"
#include "stk_mesh/base/DeviceBucket.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"

namespace stk {
namespace mesh {

using DeviceBucket = DeviceBucketT<stk::ngp::MemSpace>;

template<typename NgpMemSpace>
class DeviceMeshT : public NgpMeshBase
{
public:
  typedef NgpMemSpace ngp_mem_space;

  static_assert(Kokkos::is_memory_space_v<NgpMemSpace>);
  using MeshExecSpace     =  typename NgpMemSpace::execution_space;
  using BucketType        = DeviceBucketT<NgpMemSpace>;
  using ConnectedNodes    = typename BucketType::ConnectedNodes;
  using ConnectedEntities = typename BucketType::ConnectedEntities;
  using ConnectedOrdinals = typename BucketType::ConnectedOrdinals;
  using Permutations      = typename BucketType::Permutations;
  using MeshIndex         = FastMeshIndex;
  using ByteBuffer        = Kokkos::View<std::byte*,NgpMemSpace>;
  using HostByteBuffer    = ByteBuffer::host_mirror_type;

  KOKKOS_FUNCTION
  DeviceMeshT()
    : NgpMeshBase(),
      bulk(nullptr),
      spatial_dimension(0),
      synchronizedCount(0),
#ifndef STK_HIDE_DEPRECATED_CODE
      m_needSyncToHost(false),
#endif
      deviceMeshHostData(nullptr),
      m_meshConn()
  {}

  explicit DeviceMeshT(const stk::mesh::BulkData& b)
    : NgpMeshBase(),
      bulk(&const_cast<stk::mesh::BulkData&>(b)),
      spatial_dimension(b.mesh_meta_data().spatial_dimension()),
      synchronizedCount(0),
#ifndef STK_HIDE_DEPRECATED_CODE
      m_needSyncToHost(false),
#endif
      endRank(static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count())),
      deviceMeshHostData(nullptr),
      m_meshConn(),
      m_deviceBucketRepo(this, b.get_initial_bucket_capacity(), b.get_maximum_bucket_capacity()),
      m_deviceBufferOffsets(UnsignedViewType<NgpMemSpace>("deviceBufferOffsets", 1)),
      m_deviceMeshIndicesOffsets(UnsignedViewType<NgpMemSpace>("deviceMeshIndicesOffsets", 1)),
      m_byteBuffer("sendAndRecvData",4096),
      m_hostByteBuffer("hostSendAndRecvData", 0)
  {
    bulk->register_device_mesh();
    deviceMeshHostData = impl::get_ngp_mesh_host_data<NgpMemSpace>(*bulk);
    update();

    deviceMeshHostData->m_hostBufferOffsets = Kokkos::create_mirror_view(m_deviceBufferOffsets);;
    deviceMeshHostData->m_hostMeshIndicesOffsets = Kokkos::create_mirror_view(m_deviceMeshIndicesOffsets);
  }

  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT(const DeviceMeshT &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT(DeviceMeshT &&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT& operator=(const DeviceMeshT &) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceMeshT& operator=(DeviceMeshT &&) = default;

  KOKKOS_FUNCTION
  virtual ~DeviceMeshT() override {
#ifndef STK_HIDE_DEPRECATED_CODE
    m_needSyncToHost = false;
#endif
  }

  void update() override;

  void update_bulk_data() override;

  bool needs_update() const override {
    return synchronizedCount < bulk->synchronized_count();  // BulkData is ahead
  }

  bool needs_update_bulk_data() const override {
    return synchronizedCount > bulk->synchronized_count();  // We are ahead of BulkData
  }

  unsigned synchronized_count() const override { return synchronizedCount; }

  KOKKOS_FUNCTION
  unsigned get_spatial_dimension() const
  {
    return spatial_dimension;
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityId identifier(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()].id();
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityRank entity_rank(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()].rank();
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityKey entity_key(stk::mesh::Entity entity) const
  {
    return entityKeys[entity.local_offset()];
  }

  KOKKOS_FUNCTION
  unsigned local_id(stk::mesh::Entity entity) const
  {
    return entityLocalIds[entity.local_offset()];
  }

  KOKKOS_FUNCTION
  stk::mesh::Entity get_entity(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex& meshIndex) const
  {
    return m_deviceBucketRepo.m_buckets[rank][meshIndex.bucket_id][meshIndex.bucket_ord];
  }

  // Look up an Entity from its identifier.  This is implemented as a slow linear search and should not
  // be used for anything but testing/debugging.  If this functionality is needed in production, please
  // contact the STK Team to schedule this work.
  KOKKOS_FUNCTION
  stk::mesh::Entity linear_get_entity(stk::mesh::EntityRank rank, stk::mesh::EntityId entityId) const
  {
    for (unsigned localOffset = 0; localOffset < entityKeys.extent(0); ++localOffset) {
      auto entityKey = entityKeys[localOffset];
      if (entityKey.rank() == rank && entityKey.id() == entityId) {
        auto deviceMeshIndex = deviceMeshIndices[localOffset];
        return m_deviceBucketRepo.m_buckets[rank][deviceMeshIndex.bucket_id][deviceMeshIndex.bucket_ord];
      }
    }

    return stk::mesh::Entity();
  }

  Entity get_entity(EntityKey entityKey);

  template <typename EntityKeyView>
  EntityViewType<NgpMemSpace> get_entities(EntityKeyView entityKeyView);

  template <typename TeamMember>
  KOKKOS_FUNCTION
  Entity get_entity(TeamMember const& team, EntityKey entityKey);

  KOKKOS_FUNCTION
  ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex, stk::mesh::EntityRank connectedRank) const
  {
    auto connEntsNew = m_meshConn.get_connected_entities(get_entity(rank, entityIndex), connectedRank);
    return connEntsNew;
  }

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return get_connected_entities(rank, entityIndex, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_edges(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return get_connected_entities(rank, entityIndex, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_faces(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return get_connected_entities(rank, entityIndex, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_elements(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex) const
  {
    return get_connected_entities(rank, entityIndex, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex, stk::mesh::EntityRank connectedRank) const
  {
    auto connOrdsNew = m_meshConn.get_connected_ordinals(get_entity(rank, entityIndex), connectedRank);
    return connOrdsNew;
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_node_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_edge_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_face_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_element_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entityIndex, stk::mesh::EntityRank connectedRank) const
  {
    auto connPermsNew = m_meshConn.get_connected_permutations(get_entity(rank, entityIndex), connectedRank);
    return connPermsNew;
  }

  KOKKOS_FUNCTION
  Permutations get_node_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::NODE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_edge_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::EDGE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_face_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::FACE_RANK);
  }

  KOKKOS_FUNCTION
  Permutations get_element_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::ELEM_RANK);
  }

  KOKKOS_FUNCTION
  stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
  {
    return device_mesh_index(entity);
  }

  KOKKOS_FUNCTION
  stk::mesh::FastMeshIndex device_mesh_index(stk::mesh::Entity entity) const
  {
    return deviceMeshIndices(entity.local_offset());
  }

  stk::NgpVector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
  {
    return stk::mesh::get_bucket_ids(get_bulk_on_host(), rank, selector);
  }

  KOKKOS_FUNCTION
  EntityRank get_end_rank() const
  {
    return endRank;
  }

  KOKKOS_FUNCTION
  unsigned num_buckets(stk::mesh::EntityRank rank) const
  {
    return m_deviceBucketRepo.num_buckets(rank);
  }

  KOKKOS_FUNCTION
  const DeviceBucketT<NgpMemSpace> &get_bucket(stk::mesh::EntityRank rank, unsigned index) const
  {
    return m_deviceBucketRepo.m_buckets[rank][index];
  }

  KOKKOS_FUNCTION
  NgpCommMapIndices<NgpMemSpace> volatile_fast_shared_comm_map(stk::topology::rank_t rank, int proc,
                                                               bool includeGhosts=false) const
  {
    const size_t dataBegin = volatileFastSharedCommMapOffset[rank][proc];
    const size_t dataEnd   = includeGhosts ? volatileFastSharedCommMapOffset[rank][proc+1]
                                           : dataBegin + volatileFastSharedCommMapNumShared[rank][proc];
    NgpCommMapIndices<NgpMemSpace> buffer = Kokkos::subview(volatileFastSharedCommMap[rank],
                                                            Kokkos::pair<size_t, size_t>(dataBegin, dataEnd));
    return buffer;
  }

  template <typename... EntitiesParams>
  void batch_destroy_entities(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities)
  {
    using EntitiesMemorySpace = typename std::remove_reference<decltype(entities)>::type::memory_space;
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, EntitiesMemorySpace>::accessible,
                  "The memory space of the 'entities' View is inaccessible from the DeviceMesh execution space");

    Kokkos::View<bool*, NgpMemSpace> wasDestroyed(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "wasDestroyed"), entities.extent(0));
    batch_destroy_entities(entities, wasDestroyed);
  }

  template <typename... EntitiesParams, typename... ResultParams>
  void batch_destroy_entities(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                              const Kokkos::View<bool*, ResultParams...>& wasDestroyed)
  {
    using EntitiesMemorySpace = typename std::remove_reference<decltype(entities)>::type::memory_space;
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, EntitiesMemorySpace>::accessible,
                  "The memory space of the 'entities' View is inaccessible from the DeviceMesh execution space");
    using ResultsMemorySpace = typename std::remove_reference<decltype(wasDestroyed)>::type::memory_space;
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, ResultsMemorySpace>::accessible,
                  "The memory space of the 'wasDestroyed' View is inaccessible from the DeviceMesh execution space");
    STK_ThrowRequireMsg(wasDestroyed.extent(0) == entities.extent(0),
                        "batch_destroy_entities: 'wasDestroyed' View must be the same length as 'entities'");

    prepare_for_device_mesh_modification();

    auto destroyableEntities = impl::get_destroyable_entities(*this, entities, wasDestroyed);

    if (destroyableEntities.extent(0) == 0) {
      return;
    }

    const stk::mesh::EntityRank maxEntityRank = impl::get_max_entity_rank(*this, destroyableEntities);
    for (stk::mesh::EntityRank connectedRank = stk::topology::BEGIN_RANK;
         connectedRank < maxEntityRank; ++connectedRank) {
      m_deviceBucketRepo.batch_destroy_relations(destroyableEntities, connectedRank);
    }
    m_deviceBucketRepo.batch_destroy_entities(destroyableEntities);
    m_deviceBucketRepo.sync_from_partitions();
    increment_synchronized_count();

    impl::invalidate_destroyed_entity_maps(entityKeys, entityLocalIds, destroyableEntities, MeshExecSpace{});

    const size_t priorCount = m_entitiesDestroyedOnDevice.extent(0);
    const size_t addedCount = destroyableEntities.extent(0);
    Kokkos::resize(m_entitiesDestroyedOnDevice, priorCount + addedCount);
    Kokkos::deep_copy(Kokkos::subview(m_entitiesDestroyedOnDevice,
                                      std::pair<size_t, size_t>(priorCount, priorCount + addedCount)),
                      destroyableEntities);
  }

  template <typename... EntitiesParams>
  void batch_destroy_relations(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                               stk::mesh::EntityRank connectedRank)
  {
    using EntitiesMemorySpace = typename std::remove_reference<decltype(entities)>::type::memory_space;
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, EntitiesMemorySpace>::accessible,
                  "The memory space of the 'entities' View is inaccessible from the DeviceMesh execution space");

    prepare_for_device_mesh_modification();

    impl_batch_destroy_relations(entities, connectedRank);
  }

  template <typename... FromParams, typename... OffsetParams, typename... ToParams,
            typename... OrdinalParams, typename... PermParams>
  void batch_declare_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
                               const Kokkos::View<unsigned*, OffsetParams...>& offsets,
                               const Kokkos::View<stk::mesh::Entity*, ToParams...>& toEntities,
                               const Kokkos::View<stk::mesh::RelationIdentifier*, OrdinalParams...>& ordinals,
                               const Kokkos::View<stk::mesh::Permutation*, PermParams...>& permutations)
  {
    using FromMemorySpace = typename std::remove_reference<decltype(fromEntities)>::type::memory_space;
    using OffsetMemorySpace = typename std::remove_reference<decltype(offsets)>::type::memory_space;
    using ToMemorySpace = typename std::remove_reference<decltype(toEntities)>::type::memory_space;
    using OrdinalMemorySpace = typename std::remove_reference<decltype(ordinals)>::type::memory_space;
    using PermutationsMemorySpace = typename std::remove_reference<decltype(permutations)>::type::memory_space;

    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, FromMemorySpace>::accessible,
                  "The memory space of the 'fromEntities' View is inaccessible from the DeviceMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, OffsetMemorySpace>::accessible,
                  "The memory space of the 'offsets' View is inaccessible from the DeviceMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, ToMemorySpace>::accessible,
                  "The memory space of the 'toEntities' View is inaccessible from the DeviceMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, OrdinalMemorySpace>::accessible,
                  "The memory space of the 'ordinals' View is inaccessible from the DeviceMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, PermutationsMemorySpace>::accessible,
                  "The memory space of the 'permutations' View is inaccessible from the DeviceMesh execution space");

    STK_ThrowRequireMsg(not needs_update(),
                        "DeviceMesh cannot be modified while there are simultaneous mesh modifications in BulkData on "
                        "the host.  Please first re-acquire this DeviceMesh to automatically update it.");

    copy_all_field_data_to_device();
    sync_all_fields_to_device();

    impl_batch_declare_relations(fromEntities, offsets, toEntities, ordinals, permutations);
  }

  template <typename... FromParams, typename... OffsetParams, typename... ToParams,
            typename... OrdinalParams>
  void batch_declare_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
                               const Kokkos::View<unsigned*, OffsetParams...>& offsets,
                               const Kokkos::View<stk::mesh::Entity*, ToParams...>& toEntities,
                               const Kokkos::View<stk::mesh::RelationIdentifier*, OrdinalParams...>& ordinals)
  {
    batch_declare_relations(fromEntities, offsets, toEntities, ordinals,
                            Kokkos::View<stk::mesh::Permutation*, NgpMemSpace>());
  }

  template <typename... FromParams, typename... ToParams, typename... PermParams>
  void batch_declare_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
                               const Kokkos::View<stk::mesh::Entity**, ToParams...>& toEntities,
                               const Kokkos::View<stk::mesh::Permutation**, PermParams...>& permutations,
                               stk::mesh::EntityRank connectedRank)
  {
    Kokkos::View<unsigned*, NgpMemSpace> offsets;
    Kokkos::View<stk::mesh::Entity*, NgpMemSpace> flatToEntities;
    Kokkos::View<stk::mesh::RelationIdentifier*, NgpMemSpace> flatOrdinals;
    Kokkos::View<stk::mesh::Permutation*, NgpMemSpace> flatPermutations;
    impl_flatten_uniform_relations(fromEntities, toEntities, permutations, connectedRank,
                                   offsets, flatToEntities, flatOrdinals, flatPermutations);

    batch_declare_relations(fromEntities, offsets, flatToEntities, flatOrdinals, flatPermutations);
  }

  template <typename... FromParams, typename... ToParams>
  void batch_declare_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
                               const Kokkos::View<stk::mesh::Entity**, ToParams...>& toEntities,
                               stk::mesh::EntityRank connectedRank)
  {
    Kokkos::View<unsigned*, NgpMemSpace> offsets;
    Kokkos::View<stk::mesh::Entity*, NgpMemSpace> flatToEntities;
    Kokkos::View<stk::mesh::RelationIdentifier*, NgpMemSpace> flatOrdinals;
    Kokkos::View<stk::mesh::Permutation*, NgpMemSpace> flatPermutations;  // stays empty: no permutations supplied
    impl_flatten_uniform_relations(fromEntities, toEntities, Kokkos::View<stk::mesh::Permutation**, NgpMemSpace>(),
                                   connectedRank, offsets, flatToEntities, flatOrdinals, flatPermutations);

    batch_declare_relations(fromEntities, offsets, flatToEntities, flatOrdinals);
  }

  void clear()
  {
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<stk::topology::NUM_RANKS; rank++)
      m_deviceBucketRepo.m_buckets[rank] = BucketView();
  }

  stk::mesh::BulkData &get_bulk_on_host()
  {
    STK_ThrowRequireMsg(bulk != nullptr, "DeviceMesh::get_bulk_on_host, bulk==nullptr");
    return *bulk;
  }

  const stk::mesh::BulkData &get_bulk_on_host() const
  {
    STK_ThrowRequireMsg(bulk != nullptr, "DeviceMeshT::get_bulk_on_host, bulk==nullptr");
    return *bulk;
  }

  KOKKOS_INLINE_FUNCTION
  const impl::MeshConnectivity<stk::ngp::UVMDeviceSpace>& get_mesh_connectivity() const
  {
    return m_meshConn;
  }

  KOKKOS_INLINE_FUNCTION
  impl::MeshConnectivity<stk::ngp::UVMDeviceSpace>& get_mesh_connectivity()
  {
    return m_meshConn;
  }

  template <typename... EntityIdsParams, typename... AddPartParams, typename... EntitiesParams>
  void batch_declare_entities(stk::topology::rank_t rank,
                              const Kokkos::View<unsigned*, EntityIdsParams...>& entityIds,
                              const Kokkos::View<PartOrdinal*, AddPartParams...>& addPartOrdinals,
                              Kokkos::View<Entity*, EntitiesParams...>& requestedEntities)
  {
    using EntityIdsMemorySpace = typename std::remove_reference<decltype(entityIds)>::type::memory_space;
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, EntityIdsMemorySpace>::accessible,
                  "The memory space of the 'entities' View is inaccessible from the DeviceMesh execution space");

    impl_batch_declare_entities(rank, entityIds, addPartOrdinals, requestedEntities);
    increment_synchronized_count();
  }

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void batch_change_entity_parts(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
  {
    using EntitiesMemorySpace = typename std::remove_reference<decltype(entities)>::type::memory_space;
    using AddPartOrdinalsMemorySpace = typename std::remove_reference<decltype(addPartOrdinals)>::type::memory_space;
    using RemovePartOrdinalsMemorySpace = typename std::remove_reference<decltype(removePartOrdinals)>::type::memory_space;

    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, EntitiesMemorySpace>::accessible,
                  "The memory space of the 'entities' View is inaccessible from the DeviceMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, AddPartOrdinalsMemorySpace>::accessible,
                  "The memory space of the 'addPartOrdinals' View is inaccessible from the DeviceMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, RemovePartOrdinalsMemorySpace>::accessible,
                  "The memory space of the 'removePartOrdinals' View is inaccessible from the DeviceMesh execution space");

    prepare_for_device_mesh_modification();

    auto comm = bulk->parallel();
    stk::CommSparse commSparse(comm);
    Kokkos::View<impl::DeniedPartRemoval*, NgpMemSpace> deniedPartRemovalView;
    if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) {
      impl::require_entity_owner(*this, entities);
      impl::communicate_shared_entities(*this, commSparse, entities, addPartOrdinals, removePartOrdinals);

      int localHasRankedPartToRemove = impl::has_ranked_remove_part(get_device_bucket_repository(), removePartOrdinals) ? 1 : 0;
      int globalHasRankedPartToRemove = 0;
      stk::all_reduce_max(comm, &localHasRankedPartToRemove, &globalHasRankedPartToRemove, 1);

      if (globalHasRankedPartToRemove != 0) {
        stk::CommSparse partRemovalDenialCommSparse(comm);
        impl::communicate_denied_part_removals(*this, commSparse, partRemovalDenialCommSparse, entities, removePartOrdinals);
        deniedPartRemovalView = impl::unpack_denied_part_removals(*this, partRemovalDenialCommSparse);
      }
    }

    auto wrappedEntities = impl::wrap_entities(entities, true, false);
    impl_batch_change_entity_parts(wrappedEntities, addPartOrdinals, removePartOrdinals, deniedPartRemovalView);

    if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) {
      auto implFunc = [this](auto const& sharedWrappedEntities, auto const& sharedAddPartOrdinals, auto const& sharedRemovePartOrdinals) {
                        this->impl_batch_change_entity_parts(sharedWrappedEntities, sharedAddPartOrdinals, sharedRemovePartOrdinals, {}, false);
                      };
      impl::unpack_shared_entities_and_callback(*this, commSparse, implFunc);
    }
  }

#ifndef STK_HIDE_DEPRECATED_CODE
  STK_DEPRECATED_MSG("Use update_bulk_data() instead.")
  void sync_to_host() {
    m_needSyncToHost = false;
  }

  STK_DEPRECATED_MSG("Use need_update_bulk_data() instead.")
  bool need_sync_to_host() const override {
    return m_needSyncToHost;
  }
#endif

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void impl_batch_change_entity_parts(const Kokkos::View<impl::EntityWrapper<>*, EntitiesParams...>& entities,
                                      const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                      const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals,
                                      Kokkos::View<impl::DeniedPartRemoval*, NgpMemSpace> deniedPartRemovalView,
                                      bool expandDownward = true)
  {
    auto validEntities = impl::remove_invalid_entries_and_resize<MeshExecSpace>(entities, impl::EntityWrapper(Entity{}));

    bool hasRankedPart = impl::has_ranked_part(get_device_bucket_repository(), addPartOrdinals, removePartOrdinals);
    if (!hasRankedPart) {
      impl_batch_change_entity_parts_without_inducible_parts(validEntities, addPartOrdinals, removePartOrdinals);
    } else {
      impl_batch_change_entity_parts_with_inducible_parts(validEntities, addPartOrdinals, removePartOrdinals, deniedPartRemovalView, expandDownward);
    }
  }

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void impl_batch_change_entity_parts_with_inducible_parts(const Kokkos::View<impl::EntityWrapper<>*, EntitiesParams...>& entities,
                                                           const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                                           const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals,
                                                           Kokkos::View<impl::DeniedPartRemoval*, NgpMemSpace> deniedPartRemovalView,
                                                           bool expandDownward = true)
  {
    using PartOrdinalsViewType = typename std::remove_reference<decltype(addPartOrdinals)>::type;
    using EntityWrapperViewType = Kokkos::View<impl::EntityWrapper<>*, NgpMemSpace>;
    using NewBucketsToAddViewType = Kokkos::View<impl::NumNewBucketsToAddPerPartition*, NgpMemSpace>;
    using GrowLastBucketViewType = Kokkos::View<impl::GrowLastBucketInPartition*, NgpMemSpace>;

#ifndef STK_HIDE_DEPRECATED_CODE
    m_needSyncToHost = true;
#endif
    increment_synchronized_count();

#ifndef NDEBUG
    check_parts_are_not_internal(addPartOrdinals, removePartOrdinals);
#endif

    Kokkos::Profiling::pushRegion("construct_part_ordinal_proxy");
    EntityWrapperViewType wrappedEntities;
    if (expandDownward) {
      auto maxNumDownwardConnectedEntities = impl::get_max_num_downward_connected_entities(*this, entities);
      auto entityInterval = maxNumDownwardConnectedEntities + 1;
      auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);

      wrappedEntities = EntityWrapperViewType(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
      impl::populate_all_downward_connected_entities_and_wrap_entities(*this, entities, entityInterval, wrappedEntities);
    }
    else { // no need to populate downward connected entities if all entities for change entity parts were provided by the owner proc
      wrappedEntities = EntityWrapperViewType(Kokkos::view_alloc("wrappedEntities"), entities.extent(0));
      Kokkos::deep_copy(wrappedEntities, entities);
    }
    impl::remove_invalid_wrapped_entities_sort_unique_merge_and_resize(wrappedEntities, MeshExecSpace{});

    // determine resulting parts per entity including inducible parts
    const unsigned maxCurrentNumPartsPerEntity = impl::get_max_num_parts_per_entity(*this, wrappedEntities);
    const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

    PartOrdinalsViewType newPartOrdinalsPerEntity("newPartOrdinals", wrappedEntities.size() * maxNewNumPartsPerEntity);
    PartOrdinalsViewType sortedAddPartOrdinals = impl::get_sorted_view(addPartOrdinals);

    // create part ordinals proxy, sort and unique it (and realloc)
    // A proxy indices view to part ordinals: (rank, startPtr, length)
    using PartOrdinalsProxyViewType = Kokkos::View<impl::PartOrdinalsProxyIndices*>;
    PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc("partOrdinalsProxy", Kokkos::WithoutInitializing), wrappedEntities.size());
    impl::set_new_part_list_per_entity_with_induced_parts(*this, wrappedEntities, sortedAddPartOrdinals, removePartOrdinals,
                                                          maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy,
                                                          deniedPartRemovalView);

    PartOrdinalsProxyViewType copiedPartOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "copiedPartOrdinalsProxy"), wrappedEntities.size());
    Kokkos::deep_copy(copiedPartOrdinalsProxy, partOrdinalsProxy);

    impl::sort_and_unique_and_resize(partOrdinalsProxy, MeshExecSpace{});
    Kokkos::fence();
    Kokkos::Profiling::popRegion();

    m_deviceBucketRepo.batch_create_partitions(partOrdinalsProxy);

    using EntitySrcDestView = Kokkos::View<impl::EntitySrcDest*, NgpMemSpace>;
    EntitySrcDestView entitySrcDestView(Kokkos::view_alloc("srcDestPartitionIdPerEntity", Kokkos::WithoutInitializing), wrappedEntities.size());

    m_deviceBucketRepo.batch_get_partitions(wrappedEntities, copiedPartOrdinalsProxy, entitySrcDestView);

    NewBucketsToAddViewType numNewBucketsToAddInPartitions("NumNewBucketsToAddInPartitions", wrappedEntities.size());
    GrowLastBucketViewType growLastBucketInPartition("growLastBucketInPartition", wrappedEntities.size());
    impl::assign_dest_bucket_id_and_ordinal(*this, entitySrcDestView, numNewBucketsToAddInPartitions,
                                            growLastBucketInPartition);

    m_deviceBucketRepo.batch_grow_buckets(growLastBucketInPartition);
    m_deviceBucketRepo.batch_create_buckets(numNewBucketsToAddInPartitions);

    set_dest_bucket_ids(*this, entitySrcDestView);

    unsigned maxNumBuckets = 0;
    for (auto rank = stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank) {
      maxNumBuckets = std::max(maxNumBuckets, m_deviceBucketRepo.num_buckets(rank));
    }
    STK_ThrowAssert(maxNumBuckets != std::numeric_limits<unsigned>::max());

    m_deviceBucketRepo.batch_move_entities(entitySrcDestView);

    m_deviceBucketRepo.sync_from_partitions();

    Kokkos::fence();
  }

  void sync_all_fields_to_device() {
    const FieldVector& allFields = get_bulk_on_host().mesh_meta_data().get_fields();
    for (FieldBase* field : allFields) {
      field->sync_to_device();
    }
  }

  MeshIndexType<NgpMemSpace>& get_fast_mesh_indices() {
    return deviceMeshIndices;
  }

  const MeshIndexType<NgpMemSpace>& get_fast_mesh_indices() const {
    return deviceMeshIndices;
  }

  EntityKeyViewType<NgpMemSpace>& get_entity_keys() {
    return entityKeys;
  }

  const EntityKeyViewType<NgpMemSpace>& get_entity_keys() const {
    return entityKeys;
  }

  auto& get_ngp_parallel_sum_host_buffer_offsets() {
    return deviceMeshHostData->m_hostBufferOffsets;
  }

  auto& get_ngp_parallel_sum_host_mesh_indices_offsets() {
    return deviceMeshHostData->m_hostMeshIndicesOffsets;
  }

  auto& get_ngp_parallel_sum_device_mesh_indices_offsets() {
    return m_deviceMeshIndicesOffsets;
  }

  auto& get_ngp_parallel_sum_host_byte_buffer() {
    return m_hostByteBuffer;
  }

  auto& get_ngp_parallel_sum_device_byte_buffer() {
    return m_byteBuffer;
  }

  KOKKOS_INLINE_FUNCTION
  impl::DeviceBucketRepository<NgpMemSpace>& get_device_bucket_repository() {
    return m_deviceBucketRepo;
  }

  KOKKOS_INLINE_FUNCTION
  impl::DeviceBucketRepository<NgpMemSpace> const& get_device_bucket_repository() const {
    return m_deviceBucketRepo;
  }

  template <typename AddPartOrdinalsViewType, typename RemovePartOrdinalsViewType>
  void check_parts_are_not_internal(AddPartOrdinalsViewType const& addPartOrdinals, RemovePartOrdinalsViewType const& removePartOrdinals);

  template <typename AddPartOrdinalsViewType>
  void check_parts_are_not_internal(AddPartOrdinalsViewType const& addPartOrdinals);

  DeviceFieldDataManagerBase* get_field_data_manager(const stk::mesh::BulkData& bulk_in);

  const DeviceFieldDataManagerBase* get_field_data_manager(const stk::mesh::BulkData& bulk_in) const;


private:

  template <typename EntityIdsView, typename AddPartOrdinalsView, typename RequestedEntitiesView>
  void impl_batch_declare_entities(stk::topology::rank_t rank,
                                   const EntityIdsView& entityIds,
                                   const AddPartOrdinalsView& addPartOrdinals,
                                   RequestedEntitiesView& requestedEntities)
  {
    using PartOrdinalsViewType = typename std::remove_reference<decltype(addPartOrdinals)>::type;
    using PartOrdinalsHostViewType = typename PartOrdinalsViewType::host_mirror_type;
    using NewBucketsToAddViewType = Kokkos::View<impl::NumNewBucketsToAddPerPartition*, NgpMemSpace>;
    using GrowLastBucketViewType = Kokkos::View<impl::GrowLastBucketInPartition*, NgpMemSpace>;

#ifndef NDEBUG
    check_parts_are_not_internal(addPartOrdinals);
#endif

    const unsigned numInitialEntities = (entityKeys.extent(0) == 0) ? 0 : entityKeys.extent(0) - 1;
    const unsigned entityKeysInitialIndex = std::max(entityKeys.extent(0), 1ul);

    std::vector<PartOrdinal> allAddPartOrdinalsOnHost;
    auto addPartOrdinals_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), addPartOrdinals);
    PartVector add_parts;
    for (unsigned i = 0; i < addPartOrdinals_host.extent(0); ++i) {
      add_parts.push_back(&bulk->mesh_meta_data().get_part(addPartOrdinals_host(i)));
    }
    add_parts.push_back(bulk->mesh_meta_data().get_part("{OWNS}"));
    stk::mesh::impl::fill_add_parts_and_supersets(add_parts, allAddPartOrdinalsOnHost);

    PartOrdinalsViewType partsAndSupersets("partsAndSupersets", allAddPartOrdinalsOnHost.size());
    Kokkos::deep_copy(partsAndSupersets, PartOrdinalsHostViewType(allAddPartOrdinalsOnHost.data(), allAddPartOrdinalsOnHost.size()));

    const unsigned maxNewNumPartsPerEntity = partsAndSupersets.extent(0);
    PartOrdinalsViewType newPartOrdinalsPerEntity("newPartOrdinals", entityIds.size()*maxNewNumPartsPerEntity);
    PartOrdinalsViewType sortedAddPartOrdinals = impl::get_sorted_view(partsAndSupersets);

    // A proxy indices view to part ordinals: (rank, startPtr, length)

    if (requestedEntities.extent(0) != entityIds.extent(0)) {
      Kokkos::resize(requestedEntities, entityIds.extent(0));
    }
    Kokkos::resize(entityKeys, entityIds.extent(0) + entityKeysInitialIndex);

    Kokkos::resize(entityLocalIds, entityIds.extent(0) + entityLocalIds.extent(0));

    std::vector<size_t> requests(bulk->mesh_meta_data().entity_rank_count());
    requests[rank] = entityIds.extent(0);

    impl::fill_requested_entities_and_keys(requestedEntities, entityKeys, entityLocalIds, entityIds,
                                           rank, numInitialEntities, entityKeysInitialIndex);
    m_meshConn.update_num_entities(std::max(m_meshConn.get_num_entities(), 1ul) + entityIds.extent(0));

    // deviceMeshIndices and entityKeys are both indexed by Entity::local_offset(), so they must stay the
    // same length.  Size from entityKeys rather than growing deviceMeshIndices by entityIds.extent(0):
    // the two views are seeded from different host quantities on update() (entityKeys from
    // m_entity_keys.capacity(), deviceMeshIndices from get_size_of_entity_index_space(), i.e. its size),
    // so growing each relative to its own extent makes that capacity-vs-size skew permanent and leaves
    // the highest new local_offset out of bounds in deviceMeshIndices.
    Kokkos::resize(deviceMeshIndices, entityKeys.extent(0));

    impl::invalidate_new_entity_mesh_indices(deviceMeshIndices, requestedEntities, entityIds.extent(0));

    using PartOrdinalsProxyViewType = Kokkos::View<impl::PartOrdinalsProxyIndices*>;
    PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "partOrdinalsProxy"), entityIds.size());
    impl::set_add_part_list_per_entity(*this, rank, requestedEntities, sortedAddPartOrdinals, newPartOrdinalsPerEntity, partOrdinalsProxy);

    PartOrdinalsProxyViewType copiedPartOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "copiedPartOrdinalsProxy"), entityIds.size());
    Kokkos::deep_copy(copiedPartOrdinalsProxy, partOrdinalsProxy);

    impl::sort_and_unique_and_resize(partOrdinalsProxy, MeshExecSpace{});
    Kokkos::fence();

    m_deviceBucketRepo.batch_create_partitions(partOrdinalsProxy);

    using EntitySrcDestView = Kokkos::View<impl::EntitySrcDest*, NgpMemSpace>;
    EntitySrcDestView entitySrcDestView(Kokkos::view_alloc("srcDestPartitionIdPerEntity", Kokkos::WithoutInitializing), entityIds.size());

    m_deviceBucketRepo.batch_get_partitions(rank, requestedEntities, copiedPartOrdinalsProxy, entitySrcDestView);

    NewBucketsToAddViewType numNewBucketsToAddInPartitions("NumNewBucketsToAddInPartitions", entityIds.size());
    GrowLastBucketViewType growLastBucketInPartition("growLastBucketInPartition", entityIds.size());
    impl::assign_dest_bucket_id_and_ordinal(*this, entitySrcDestView, numNewBucketsToAddInPartitions,
                                            growLastBucketInPartition);

    {
      const MetaData& meta = get_bulk_on_host().mesh_meta_data();
      const FieldVector& fields = meta.get_fields();
      for(FieldBase* field : fields) {
        field_datatype_execute(*field,
          [&]<typename T>(const stk::mesh::FieldBase& fieldBase) {
            fieldBase.data<T, ReadWrite, stk::ngp::DeviceSpace>();
          }
        );
      }
    }

    m_deviceBucketRepo.batch_grow_buckets(growLastBucketInPartition);
    m_deviceBucketRepo.batch_create_buckets(numNewBucketsToAddInPartitions);

    set_dest_bucket_ids(*this, entitySrcDestView);

    unsigned maxNumBuckets = 0;
    for (auto current_rank = stk::topology::BEGIN_RANK; current_rank < stk::topology::END_RANK; ++current_rank) {
      maxNumBuckets = std::max(maxNumBuckets, m_deviceBucketRepo.num_buckets(current_rank));
    }
    STK_ThrowAssert(maxNumBuckets != std::numeric_limits<unsigned>::max());

    m_deviceBucketRepo.batch_declare_entities(entitySrcDestView);

    m_deviceBucketRepo.sync_from_partitions();

    for (stk::mesh::EntityRank rankIdx = stk::topology::NODE_RANK; rankIdx < endRank; ++rankIdx) {
      for(unsigned i=0; i<m_deviceBucketRepo.num_buckets(rankIdx); ++i) {
        m_deviceBucketRepo.get_bucket(rankIdx,i)->set_mesh_connectivity(m_meshConn);
      }
    }

    Kokkos::fence();

    synchronizedCount++;
  }

  template <typename... EntitiesParams>
  void impl_batch_destroy_relations(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                                    stk::mesh::EntityRank connectedRank)
  {
    m_deviceBucketRepo.batch_destroy_relations(entities, connectedRank);
    m_deviceBucketRepo.sync_from_partitions();
    increment_synchronized_count();
  }

  template <typename... FromParams, typename... OffsetParams, typename... ToParams,
            typename... OrdinalParams, typename... PermParams>
  void impl_batch_declare_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
                                    const Kokkos::View<unsigned*, OffsetParams...>& offsets,
                                    const Kokkos::View<stk::mesh::Entity*, ToParams...>& toEntities,
                                    const Kokkos::View<stk::mesh::RelationIdentifier*, OrdinalParams...>& ordinals,
                                    const Kokkos::View<stk::mesh::Permutation*, PermParams...>& permutations)
  {
    m_deviceBucketRepo.batch_declare_relations(fromEntities, offsets, toEntities, ordinals, permutations);
    m_deviceBucketRepo.sync_from_partitions();
    increment_synchronized_count();
  }

  template <typename... FromParams, typename... ToParams, typename... PermParams>
  void impl_flatten_uniform_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
                                      const Kokkos::View<stk::mesh::Entity**, ToParams...>& toEntities,
                                      const Kokkos::View<stk::mesh::Permutation**, PermParams...>& permutations,
                                      [[maybe_unused]] stk::mesh::EntityRank connectedRank,
                                      Kokkos::View<unsigned*, NgpMemSpace>& offsets,
                                      Kokkos::View<stk::mesh::Entity*, NgpMemSpace>& flatToEntities,
                                      Kokkos::View<stk::mesh::RelationIdentifier*, NgpMemSpace>& flatOrdinals,
                                      Kokkos::View<stk::mesh::Permutation*, NgpMemSpace>& flatPermutations)
  {
    const size_t numFrom = fromEntities.extent(0);
    const size_t numCols = toEntities.extent(1);
    const bool hasPermutations = permutations.extent(0) > 0;

    auto toEntitiesHost   = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, toEntities);
    auto permutationsHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, permutations);

#ifndef NDEBUG
    if (numFrom > 0) {
      // Only the debug validation inspects the from-entities on host; skip its mirror in release.
      auto fromEntitiesHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, fromEntities);
      STK_ThrowRequire(toEntities.extent(0) == numFrom);   // one row of connectivity per from-entity
      STK_ThrowRequire(bulk->is_valid(fromEntitiesHost(0)));
      const stk::topology topo = bulk->bucket(fromEntitiesHost(0)).topology();
      const unsigned expectedCols = topo.num_sub_topology(connectedRank);
      STK_ThrowRequire(expectedCols > 0);          // connectedRank must be a valid downward sub-rank
      STK_ThrowRequire(numCols == expectedCols);   // full complement required
      for (size_t i = 0; i < numFrom; ++i) {
        STK_ThrowRequire(bulk->is_valid(fromEntitiesHost(i)));
        STK_ThrowRequire(bulk->bucket(fromEntitiesHost(i)).topology() == topo);  // single uniform topology
      }
      if (hasPermutations) {   // permutations grid must mirror the toEntities grid
        STK_ThrowRequire(permutations.extent(0) == numFrom);
        STK_ThrowRequire(permutations.extent(1) == numCols);
      }
    }
#endif

    const size_t numRelations = numFrom * numCols;
    Kokkos::resize(offsets, numFrom > 0 ? numFrom + 1 : 0);
    Kokkos::resize(flatToEntities, numRelations);
    Kokkos::resize(flatOrdinals, numRelations);
    if (hasPermutations) {
      Kokkos::resize(flatPermutations, numRelations);
    }

    auto offsetsHost = Kokkos::create_mirror_view(offsets);
    auto flatToEntitiesHost  = Kokkos::create_mirror_view(flatToEntities);
    auto flatOrdinalsHost = Kokkos::create_mirror_view(flatOrdinals);
    auto flatPermutationsHost = Kokkos::create_mirror_view(flatPermutations);

    for (size_t i = 0; i < numFrom; ++i) {
      const size_t base = i * numCols;
      offsetsHost(i) = static_cast<unsigned>(base);
      for (size_t j = 0; j < numCols; ++j) {
        flatToEntitiesHost(base + j) = toEntitiesHost(i, j);
        flatOrdinalsHost(base + j) = static_cast<stk::mesh::RelationIdentifier>(j);
        if (hasPermutations) {
          flatPermutationsHost(base + j) = permutationsHost(i, j);
        }
      }
    }
    if (numFrom > 0) {
      offsetsHost(numFrom) = static_cast<unsigned>(numRelations);
    }

    Kokkos::deep_copy(offsets, offsetsHost);
    Kokkos::deep_copy(flatToEntities, flatToEntitiesHost);
    Kokkos::deep_copy(flatOrdinals, flatOrdinalsHost);
    if (hasPermutations) {
      Kokkos::deep_copy(flatPermutations, flatPermutationsHost);
    }
  }

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void impl_batch_change_entity_parts_without_inducible_parts(const Kokkos::View<impl::EntityWrapper<>*, EntitiesParams...>& entities,
                                                              const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                                              const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
  {
    using PartOrdinalsViewType = typename std::remove_reference<decltype(addPartOrdinals)>::type;
    using NewBucketsToAddViewType = Kokkos::View<impl::NumNewBucketsToAddPerPartition*, NgpMemSpace>;
    using GrowLastBucketViewType = Kokkos::View<impl::GrowLastBucketInPartition*, NgpMemSpace>;

#ifndef STK_HIDE_DEPRECATED_CODE
    m_needSyncToHost = true;
#endif
    increment_synchronized_count();

#ifndef NDEBUG
    check_parts_are_not_internal(addPartOrdinals, removePartOrdinals);
#endif

    const unsigned maxCurrentNumPartsPerEntity = impl::get_max_num_parts_per_entity(*this, entities);
    const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

    PartOrdinalsViewType newPartOrdinalsPerEntity("newPartOrdinals", entities.size()*maxNewNumPartsPerEntity);
    PartOrdinalsViewType sortedAddPartOrdinals = impl::get_sorted_view(addPartOrdinals);

    // A proxy indices view to part ordinals: (rank, startPtr, length)
    using PartOrdinalsProxyViewType = Kokkos::View<impl::PartOrdinalsProxyIndices*>;
    PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "partOrdinalsProxy"), entities.size());
    impl::set_new_part_list_per_entity(*this, entities, sortedAddPartOrdinals, removePartOrdinals,
                                       maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);

    PartOrdinalsProxyViewType copiedPartOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "copiedPartOrdinalsProxy"), entities.size());
    Kokkos::deep_copy(copiedPartOrdinalsProxy, partOrdinalsProxy);

    impl::sort_and_unique_and_resize(partOrdinalsProxy, MeshExecSpace{});
    Kokkos::fence();

    m_deviceBucketRepo.batch_create_partitions(partOrdinalsProxy);

    using EntitySrcDestView = Kokkos::View<impl::EntitySrcDest*, NgpMemSpace>;
    EntitySrcDestView entitySrcDestView(Kokkos::view_alloc("srcDestPartitionIdPerEntity", Kokkos::WithoutInitializing), entities.size());

    m_deviceBucketRepo.batch_get_partitions(entities, copiedPartOrdinalsProxy, entitySrcDestView);

    NewBucketsToAddViewType numNewBucketsToAddInPartitions("NumNewBucketsToAddInPartitions", entities.size());
    GrowLastBucketViewType growLastBucketInPartition("growLastBucketInPartition", entities.size());
    impl::assign_dest_bucket_id_and_ordinal(*this, entitySrcDestView, numNewBucketsToAddInPartitions,
                                            growLastBucketInPartition);

    m_deviceBucketRepo.batch_grow_buckets(growLastBucketInPartition);
    m_deviceBucketRepo.batch_create_buckets(numNewBucketsToAddInPartitions);

    set_dest_bucket_ids(*this, entitySrcDestView);

    unsigned maxNumBuckets = 0;
    for (auto rank = stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank) {
      maxNumBuckets = std::max(maxNumBuckets, m_deviceBucketRepo.num_buckets(rank));
    }
    STK_ThrowAssert(maxNumBuckets != std::numeric_limits<unsigned>::max());

    m_deviceBucketRepo.batch_move_entities(entitySrcDestView);

    m_deviceBucketRepo.sync_from_partitions();

    Kokkos::fence();
  }

  bool fill_buckets(const stk::mesh::BulkData& bulk_in);

  bool update_mesh_connectivity(const stk::mesh::BulkData& bulk_in);

  void copy_entity_keys_to_device();

  void copy_entity_local_ids_to_device();

  void copy_volatile_fast_shared_comm_map_to_device();

  void update_field_data_manager();

  void update_field_metadata_host_pointers();

  void increment_synchronized_count() { ++synchronizedCount; }

  void prepare_for_device_mesh_modification()
  {
    STK_ThrowRequireMsg(not needs_update(),
                        "DeviceMesh cannot be modified while there are simultaneous mesh modifications in BulkData on "
                        "the host.  Please first re-acquire this DeviceMesh to automatically update it.");
    copy_all_field_data_to_device();
    sync_all_fields_to_device();
  }

  void copy_all_field_data_to_device();

  using BucketView = Kokkos::View<DeviceBucketT<NgpMemSpace>*, stk::ngp::UVMMemSpace>;
  stk::mesh::BulkData* bulk;
  unsigned spatial_dimension;
  unsigned synchronizedCount;

  Kokkos::View<stk::mesh::Entity*, NgpMemSpace> m_entitiesDestroyedOnDevice;

#ifndef STK_HIDE_DEPRECATED_CODE
  bool m_needSyncToHost;
#endif

  stk::mesh::EntityRank endRank;
  impl::NgpMeshHostData<NgpMemSpace>* deviceMeshHostData;
  impl::MeshConnectivity<stk::ngp::UVMDeviceSpace> m_meshConn;

  EntityKeyViewType<NgpMemSpace> entityKeys;
  UnsignedViewType<NgpMemSpace> entityLocalIds;

  impl::DeviceBucketRepository<NgpMemSpace> m_deviceBucketRepo;
  HostMeshIndexType<NgpMemSpace> hostMeshIndices;
  MeshIndexType<NgpMemSpace> deviceMeshIndices;

  UnsignedViewType<NgpMemSpace> volatileFastSharedCommMapOffset[stk::topology::NUM_RANKS];
  UnsignedViewType<NgpMemSpace> volatileFastSharedCommMapNumShared[stk::topology::NUM_RANKS];
  FastSharedCommMapViewType<NgpMemSpace> volatileFastSharedCommMap[stk::topology::NUM_RANKS];

  UnsignedViewType<NgpMemSpace> m_deviceBufferOffsets;
  UnsignedViewType<NgpMemSpace> m_deviceMeshIndicesOffsets;
  ByteBuffer m_byteBuffer;
  HostByteBuffer m_hostByteBuffer;
};

using DeviceMesh = DeviceMeshT<stk::ngp::MemSpace>;

constexpr double RESIZE_FACTOR = 0.05;

template <typename DEVICE_VIEW, typename HOST_VIEW>
inline void reallocate_views(DEVICE_VIEW & deviceView, HOST_VIEW & hostView, size_t requiredSize, double resizeFactor = 0.0)
{
  const size_t currentSize = deviceView.extent(0);
  const size_t shrinkThreshold = currentSize - static_cast<size_t>(2*resizeFactor*currentSize);
  const bool needGrowth = (requiredSize > currentSize);
  const bool needShrink = (requiredSize < shrinkThreshold);

  if (needGrowth || needShrink) {
    const size_t newSize = requiredSize + static_cast<size_t>(resizeFactor*requiredSize);
    deviceView = DEVICE_VIEW(Kokkos::view_alloc(Kokkos::WithoutInitializing, deviceView.label()), newSize);
    hostView = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, deviceView);
  }
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::update()
{
  STK_ThrowRequireMsg(not bulk->in_modifiable_state(),
                      "Cannot update DeviceMesh while the host BulkData is being modified.");

  if (not needs_update()) {
    return;
  }

  require_ngp_mesh_rank_limit(bulk->mesh_meta_data());

  Kokkos::Profiling::pushRegion("DeviceMeshT::update_mesh");
  const bool anyChanges = fill_buckets(*bulk);

  if (anyChanges) {
    Kokkos::Profiling::pushRegion("anyChanges stuff");

    Kokkos::Profiling::pushRegion("entity-keys");
    copy_entity_keys_to_device();
    copy_entity_local_ids_to_device();
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("mesh-indices");
    deviceMeshIndices = bulk->get_updated_fast_mesh_indices<NgpMemSpace>();
    Kokkos::Profiling::popRegion();

    update_field_data_manager();

    Kokkos::Profiling::popRegion();
  }

  Kokkos::Profiling::pushRegion("volatile-fast-shared-comm-map");
  copy_volatile_fast_shared_comm_map_to_device();
  Kokkos::Profiling::popRegion();

  synchronizedCount = bulk->synchronized_count();
  Kokkos::Profiling::popRegion();
}

template <typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::update_bulk_data()
{
#ifndef STK_HIDE_DEPRECATED_CODE
  m_needSyncToHost = false;
#endif

  STK_ThrowRequireMsg(not bulk->in_modifiable_state(),
                      "Cannot update BulkData from DeviceMesh while it is already being modified.");

  if (not needs_update_bulk_data()) {
    return;
  }

  // The mesh modifications applied below will automatically increment the fieldMetaDataModCount value,
  // putting it ahead of the mod state on device and triggering another update back to device.  To
  // avoid this infinite loop, store off the initial global mod counts (which match device), update the mesh,
  // and reset the values to match device since the mesh actually does match now.
  const stk::mesh::FieldVector& allFields = bulk->mesh_meta_data().get_fields();
  std::vector<unsigned> cachedFieldMetaDataModCounts(allFields.size());
  for (unsigned i = 0; i < allFields.size(); ++i) {
    cachedFieldMetaDataModCounts[i] = allFields[i]->data_bytes<const std::byte>().field_meta_data_mod_count();
  }

  bulk->modification_begin_for_sync_to_host("sync DeviceMesh to host");

  m_deviceBucketRepo.sync_to_host(bulk->m_bucket_repository);

  while(bulk->get_size_of_entity_index_space() < entityKeys.extent(0)) {
    bulk->generate_new_entity(0u);
  }

  // Order dependency: this block MUST run before the entityKeys -> bulk->m_entity_keys deep_copy below.
  // record_entity_deletion reads entity_key(entity) (from the host m_entity_keys) to unregister the entity
  // from the key mapping, but batch_destroy_entities has already invalidated these entities' slots in the
  // device entityKeys view.  If the deep_copy ran first it would overwrite the host keys with those
  // invalidated values, so record_entity_deletion would read an empty key and fail to unregister.
  if (m_entitiesDestroyedOnDevice.extent(0) > 0) {
    auto hostDestroyedEntities = Kokkos::create_mirror_view(m_entitiesDestroyedOnDevice);
    Kokkos::deep_copy(hostDestroyedEntities, m_entitiesDestroyedOnDevice);
    for (size_t i = 0; i < hostDestroyedEntities.extent(0); ++i) {
      const stk::mesh::Entity destroyedEntity = hostDestroyedEntities(i);
      if (bulk->is_valid(destroyedEntity)) {
        bulk->record_entity_deletion(destroyedEntity, false);
      }
    }
    Kokkos::resize(m_entitiesDestroyedOnDevice, 0);
  }

  // Copies the (post-destroy, invalidated) device keys onto the host; see the ordering note above.
  Kokkos::deep_copy(Kokkos::subview(bulk->m_entity_keys.get_view(), std::pair<int, int>(1, entityKeys.extent(0))),
                    Kokkos::subview(entityKeys, std::pair<int, int>(1, entityKeys.extent(0))));

  for (int i=0; i < stk::topology::NUM_RANKS; ++i)
  {
    Kokkos::resize(Kokkos::WithoutInitializing, deviceMeshHostData->hostVolatileFastSharedCommMap[i],
                   volatileFastSharedCommMap[i].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, deviceMeshHostData->hostVolatileFastSharedCommMapOffset[i],
                   volatileFastSharedCommMapOffset[i].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, deviceMeshHostData->hostVolatileFastSharedCommMapNumShared[i],
                   volatileFastSharedCommMapNumShared[i].extent(0));

    Kokkos::deep_copy(deviceMeshHostData->hostVolatileFastSharedCommMap[i], volatileFastSharedCommMap[i]);
    Kokkos::deep_copy(deviceMeshHostData->hostVolatileFastSharedCommMapOffset[i], volatileFastSharedCommMapOffset[i]);
    Kokkos::deep_copy(deviceMeshHostData->hostVolatileFastSharedCommMapNumShared[i],
                      volatileFastSharedCommMapNumShared[i]);
  }

  bulk->m_meshModification.set_sync_count(synchronizedCount);
  bulk->modification_end_for_sync_to_host();

  for (EntityRank rank=stk::topology::BEGIN_RANK; rank < bulk->mesh_meta_data().entity_rank_count(); ++rank)
  {
    for (Bucket* bucket : bulk->buckets(rank))
    {
      for (unsigned i=0; i < bucket->size(); ++i)
      {
        bulk->set_mesh_index((*bucket)[i], bucket, i);
      }
    }
  }

  update_field_metadata_host_pointers();

  for (unsigned i = 0; i < allFields.size(); ++i) {
    auto& fieldDataBytes = allFields[i]->data_bytes<const std::byte>();
    fieldDataBytes.set_field_meta_data_mod_count(cachedFieldMetaDataModCounts[i]);
    fieldDataBytes.set_up_to_date();
  }
}

template<typename NgpMemSpace>
bool DeviceMeshT<NgpMemSpace>::fill_buckets(const stk::mesh::BulkData& bulk_in)
{
  bool anyBucketChanges = false;
  update_mesh_connectivity(bulk_in);

  Kokkos::Profiling::pushRegion("fill_buckets");
  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
    auto& hostBuckets = bulk_in.buckets(rank);
    if (static_cast<unsigned>(rank) < bulk_in.m_bucket_repository.m_partitions.size()) {
      auto& hostPartitions = bulk_in.m_bucket_repository.m_partitions[rank];
      m_deviceBucketRepo.copy_buckets_and_partitions_from_host(rank, hostBuckets, hostPartitions, anyBucketChanges);
    }
  }

  m_meshConn.set_update_count(bulk_in.synchronized_count());

  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
    for(unsigned i=0; i<m_deviceBucketRepo.num_buckets(rank); ++i) {
      m_deviceBucketRepo.get_bucket(rank,i)->set_mesh_connectivity(m_meshConn);
    }
  }

  Kokkos::Profiling::popRegion();

  return anyBucketChanges;
}

template<typename NgpMemSpace>
bool DeviceMeshT<NgpMemSpace>::update_mesh_connectivity(const stk::mesh::BulkData& bulk_in)
{
  Kokkos::Profiling::pushRegion("update_mesh_connectivity");
  bool returnValue = false;
  const size_t numBulkEntities = bulk_in.get_size_of_entity_index_space();
  const size_t numMeshConnEntities = m_meshConn.get_num_entities();
  const size_t diff = static_cast<size_t>(std::abs(static_cast<int64_t>(numBulkEntities - numMeshConnEntities)));
  const size_t percentIncrease = numMeshConnEntities>0 ?
    static_cast<size_t>(((1.0*diff)/numMeshConnEntities)*100) : 100;
  bool needToReplaceConnectivity = percentIncrease > 20;

  if (!needToReplaceConnectivity) {
    const size_t numMeshModsSinceLastUpdate = bulk_in.synchronized_count() - m_meshConn.get_update_count();
    auto [allocBlocksNeeded, modifiedEntities] = (numMeshModsSinceLastUpdate==1) ?
      impl::dealloc_and_get_needed_allocs_using_entity_states(bulk_in, m_meshConn) :
      impl::dealloc_and_get_needed_allocs_full_mesh(bulk_in, m_meshConn);

    auto allocBlocksAvail = impl::get_available_blocks(m_meshConn.get_memory_pool());

    bool canUpdateConnectivity = impl::blocks_are_available(allocBlocksNeeded, allocBlocksAvail);
    if (canUpdateConnectivity) {
      impl::update_mesh_connectivity(bulk_in, m_meshConn, modifiedEntities);
      returnValue = true;
    }
    else {
      needToReplaceConnectivity = true;
    }
  }

  if (needToReplaceConnectivity) {
    impl::fill_mesh_connectivity(bulk_in, m_meshConn);
    returnValue = true;
  }

  Kokkos::Profiling::popRegion();
  return returnValue;
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_entity_keys_to_device()
{
  const size_t hostCapacity = get_bulk_on_host().m_entity_keys.capacity();
  const size_t hostSize = get_bulk_on_host().m_entity_keys.size();

  if (hostCapacity != entityKeys.extent(0)) {
    Kokkos::resize(Kokkos::WithoutInitializing, entityKeys, hostCapacity);
  }

  Kokkos::deep_copy(Kokkos::subview(entityKeys, std::pair<size_t, size_t>(0, hostSize)),
                    Kokkos::subview(get_bulk_on_host().m_entity_keys.get_view(),
                                    std::pair<size_t, size_t>(0, hostSize)));

  if (hostSize < hostCapacity) {
    Kokkos::deep_copy(Kokkos::subview(entityKeys, std::pair<size_t, size_t>(hostSize, hostCapacity)),
                      stk::mesh::EntityKey());
  }
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_entity_local_ids_to_device()
{
  if (get_bulk_on_host().m_local_ids.capacity() != entityLocalIds.extent(0)) {
    Kokkos::resize(Kokkos::WithoutInitializing, entityLocalIds, get_bulk_on_host().m_local_ids.capacity());
  }

  Kokkos::deep_copy(entityLocalIds, get_bulk_on_host().m_local_ids.get_view());
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_volatile_fast_shared_comm_map_to_device()
{
  bulk->volatile_fast_shared_comm_map<NgpMemSpace>(stk::topology::NODE_RANK, 0);
  auto& hostVolatileFastSharedCommMapOffset = deviceMeshHostData->hostVolatileFastSharedCommMapOffset;
  auto& hostVolatileFastSharedCommMapNumShared = deviceMeshHostData->hostVolatileFastSharedCommMapNumShared;
  auto& hostVolatileFastSharedCommMap = deviceMeshHostData->hostVolatileFastSharedCommMap;

  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEM_RANK; ++rank)
  {
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMapNumShared[rank], hostVolatileFastSharedCommMapNumShared[rank].extent(0));
    Kokkos::resize(Kokkos::WithoutInitializing, volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank].extent(0));
    Kokkos::deep_copy(volatileFastSharedCommMapOffset[rank], hostVolatileFastSharedCommMapOffset[rank]);
    Kokkos::deep_copy(volatileFastSharedCommMapNumShared[rank], hostVolatileFastSharedCommMapNumShared[rank]);
    Kokkos::deep_copy(volatileFastSharedCommMap[rank], hostVolatileFastSharedCommMap[rank]);
  }
}

template<typename NgpMemSpace>
DeviceFieldDataManagerBase* DeviceMeshT<NgpMemSpace>::get_field_data_manager(const stk::mesh::BulkData& bulk_in)
{
  DeviceFieldDataManagerBase* deviceFieldDataManagerBase = nullptr;

  if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::DeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::DeviceSpace>();
  }
  else if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::UVMDeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::UVMDeviceSpace>();
  }
  else if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::HostPinnedDeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::HostPinnedDeviceSpace>();
  }
  else {
    STK_ThrowErrorMsg("Requested a DeviceFieldDataManager from a DeviceMesh with an unsupported MemorySpace: " <<
                      typeid(NgpMemSpace).name());
  }

  return deviceFieldDataManagerBase;
}

template<typename NgpMemSpace>
const DeviceFieldDataManagerBase* DeviceMeshT<NgpMemSpace>::get_field_data_manager(const stk::mesh::BulkData& bulk_in) const
{
  DeviceFieldDataManagerBase* deviceFieldDataManagerBase = nullptr;

  if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::DeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::DeviceSpace>();
  }
  else if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::UVMDeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::UVMDeviceSpace>();
  }
  else if constexpr (std::is_same_v<NgpMemSpace, stk::ngp::HostPinnedDeviceSpace::mem_space>) {
    deviceFieldDataManagerBase = bulk_in.get_device_field_data_manager<stk::ngp::HostPinnedDeviceSpace>();
  }
  else {
    STK_ThrowErrorMsg("Requested a DeviceFieldDataManager from a DeviceMesh with an unsupported MemorySpace: " <<
                      typeid(NgpMemSpace).name());
  }

  return deviceFieldDataManagerBase;
}

template<typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::update_field_data_manager()
{
  DeviceFieldDataManagerBase* deviceFieldDataManagerBase = get_field_data_manager(*bulk);
  STK_ThrowRequire(deviceFieldDataManagerBase != nullptr);
  deviceFieldDataManagerBase->update_all_bucket_allocations();
}

template<typename NgpMemSpace>
template <typename AddPartOrdinalsViewType, typename RemovePartOrdinalsViewType>
void DeviceMeshT<NgpMemSpace>::check_parts_are_not_internal(AddPartOrdinalsViewType const& addPartOrdinals,
                                                            RemovePartOrdinalsViewType const& removePartOrdinals)
{
  auto& meta = bulk->mesh_meta_data();

  auto addPartOrdinalsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, addPartOrdinals);
  for (unsigned i = 0; i < addPartOrdinalsHost.extent(0); ++i) {
    if (impl::is_internal_part(meta.get_part(addPartOrdinalsHost(i)))) {
      Kokkos::abort("Cannot add an internal part.\n");
    }
  }

  auto removePartOrdinalsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, removePartOrdinals);
  for (unsigned i = 0; i < removePartOrdinalsHost.extent(0); ++i) {
    if (impl::is_internal_part(meta.get_part(removePartOrdinalsHost(i)))) {
      Kokkos::abort("Cannot remove an internal part.\n");
    }
  }
}

template<typename NgpMemSpace>
template <typename AddPartOrdinalsViewType>
void DeviceMeshT<NgpMemSpace>::check_parts_are_not_internal(AddPartOrdinalsViewType const& addPartOrdinals)
{
  auto& meta = bulk->mesh_meta_data();

  auto addPartOrdinalsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, addPartOrdinals);
  for (unsigned i = 0; i < addPartOrdinalsHost.extent(0); ++i) {
    if (impl::is_internal_part(meta.get_part(addPartOrdinalsHost(i)))) {
      Kokkos::abort("Cannot add an internal part.\n");
    }
  }
}

template <typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::update_field_metadata_host_pointers()
{
  {
    for (FieldBase* field : bulk->mesh_meta_data().get_fields())
    {
      get_field_data_manager(*bulk)->update_host_bucket_pointers(field->mesh_meta_data_ordinal());
      field->modify_on_device();
    }
  }
}

template <typename NgpMemSpace>
void DeviceMeshT<NgpMemSpace>::copy_all_field_data_to_device()
{
  const FieldVector& allFields = get_bulk_on_host().mesh_meta_data().get_fields();

  for (auto field : allFields) {
    if (!field->has_device_data()) {
      field_datatype_execute(*field,
        [&]<typename T>(const stk::mesh::FieldBase& fieldBase) {
          fieldBase.data<T,stk::mesh::ReadOnly,stk::ngp::DeviceSpace>();
        }
      );
    }
  }
}

template <typename NgpMemSpace>
Entity DeviceMeshT<NgpMemSpace>::get_entity(EntityKey entityKey)
{
  STK_ThrowRequire(entityKey.is_valid());

  using OffsetType = stk::mesh::Entity::entity_value_type;
  OffsetType foundOffset = 0;
  // Clamp to the offsets that both views actually cover, matching get_entities() and the team-level
  // get_entity() overload below.  entityKeys is sized from the host m_entity_keys capacity, which can
  // exceed the entity index space.
  const size_t numEntityOffsets = Kokkos::min(entityKeys.extent(0), deviceMeshIndices.extent(0));
  Kokkos::parallel_reduce(numEntityOffsets,
    KOKKOS_CLASS_LAMBDA(const int i, OffsetType& update) {
      if (entityKeys(i) == entityKey) {
        update = static_cast<OffsetType>(i);
      }
    }, Kokkos::Max<OffsetType>(foundOffset)
  );

  return Entity{foundOffset};
}

template <typename NgpMemSpace>
template <typename EntityKeyView>
EntityViewType<NgpMemSpace> DeviceMeshT<NgpMemSpace>::get_entities(EntityKeyView entityKeyView)
{
  using TeamPolicy = Kokkos::TeamPolicy<MeshExecSpace>;
  using TeamMember = typename TeamPolicy::member_type;

  auto numEntityKeys = entityKeyView.extent(0);
  Kokkos::View<Entity*> tempEntities("tempEntities", numEntityKeys);

  unsigned numValidEntities = 0;

  auto& entityKeysInMesh = entityKeys;
  const size_t numMeshOffsets = Kokkos::min(entityKeysInMesh.extent(0), deviceMeshIndices.extent(0));
  Kokkos::parallel_reduce(TeamPolicy(numMeshOffsets, Kokkos::AUTO),
    KOKKOS_CLASS_LAMBDA(TeamMember const& team, unsigned& update) {
      auto i = team.league_rank();
      auto entityKeyInMesh = entityKeysInMesh(i);
      auto found = false;

      if (!entityKeyInMesh.is_valid()) { return; }

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, entityKeyView.extent(0)),
        [&](const int j) {
          auto localKey = entityKeyView(j);

          if (found) { return; }

          if (localKey == entityKeyInMesh) {
            STK_NGP_ThrowAssert(!found);
            tempEntities(j) = Entity{static_cast<unsigned>(i)};
            found = true;
          }
        });
      if (found) ++update;
    }, Kokkos::Sum<unsigned>(numValidEntities)
  );
  Kokkos::fence();

  if (numValidEntities == numEntityKeys) {
    return tempEntities;
  }

  Kokkos::View<Entity*> compactEntities("entities", numValidEntities);
  Kokkos::parallel_scan(numEntityKeys,
    KOKKOS_LAMBDA(const int i, unsigned& update, bool isFinal) {
      bool isValid = tempEntities(i).is_local_offset_valid();
      unsigned pos = update;

      if (isValid) { ++update; }
      if (isFinal && isValid) { compactEntities(pos) = tempEntities(i); }
    }
  );
  Kokkos::fence();

  return compactEntities;
}

template <typename NgpMemSpace>
template <typename TeamMember>
KOKKOS_FUNCTION
Entity DeviceMeshT<NgpMemSpace>::get_entity(TeamMember const& team, EntityKey entityKey)
{
  using OffsetType = Entity::entity_value_type;

  auto& entityKeyInMesh = entityKeys;
  OffsetType foundOffset = 0;
  auto numEntityOffsets = Kokkos::min(entityKeyInMesh.extent(0), deviceMeshIndices.extent(0));

  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, numEntityOffsets),
    [&](const unsigned i, OffsetType& update) {
      update = (entityKeyInMesh(i) == entityKey) ? static_cast<OffsetType>(i) : OffsetType{};
    }, Kokkos::Max<OffsetType>(foundOffset)
  );

  return Entity{foundOffset};
}

}
}

#endif

