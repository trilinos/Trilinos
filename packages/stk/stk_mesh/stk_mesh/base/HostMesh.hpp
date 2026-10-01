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

#ifndef STK_MESH_HOSTMESH_HPP
#define STK_MESH_HOSTMESH_HPP

#include <stk_util/stk_config.h>
#include <stk_util/util/StridedArray.hpp>
#include "stk_mesh/base/NgpMeshBase.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/baseImpl/BucketRepository.hpp"
#include "stk_mesh/base/DestroyRelations.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_topology/topology.hpp"
#include <Kokkos_Core.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/NgpUtils.hpp>
#include <stk_util/util/StkNgpVector.hpp>

namespace stk {
namespace mesh {

template<typename NgpMemSpace>
class HostMeshT : public NgpMeshBase
{
public:
  typedef NgpMemSpace ngp_mem_space;

  static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, NgpMemSpace>::accessible);
  static_assert(Kokkos::is_memory_space_v<NgpMemSpace>);
  using MeshExecSpace     = typename NgpMemSpace::execution_space;
  using MeshIndex         = FastMeshIndex;
  using BucketType        = stk::mesh::Bucket;
  using ConnectedNodes    = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedOrdinals = util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
  using Permutations      = util::StridedArray<const stk::mesh::Permutation>;

  KOKKOS_FUNCTION
  HostMeshT()
    : NgpMeshBase(),
      bulk(nullptr)
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
  }

  KOKKOS_FUNCTION
  HostMeshT(const stk::mesh::BulkData& b)
    : NgpMeshBase(),
      bulk(&const_cast<stk::mesh::BulkData&>(b))
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(m_syncCountWhenUpdated = bulk->synchronized_count(););
    KOKKOS_IF_ON_HOST(require_ngp_mesh_rank_limit(bulk->mesh_meta_data()););
  }

  KOKKOS_FUNCTION virtual ~HostMeshT() override {}

  HostMeshT(const HostMeshT &) = default;
  HostMeshT(HostMeshT &&) = default;
  HostMeshT& operator=(const HostMeshT&) = default;
  HostMeshT& operator=(HostMeshT&&) = default;

  void update() override {
    m_syncCountWhenUpdated = bulk->synchronized_count();
  }

  void update_bulk_data() override {}

  bool needs_update() const override {
    return m_syncCountWhenUpdated != bulk->synchronized_count();
  }

  bool needs_update_bulk_data() const override { return false; }

  unsigned synchronized_count() const override { return m_syncCountWhenUpdated; }

  KOKKOS_FUNCTION
  unsigned get_spatial_dimension() const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return bulk->mesh_meta_data().spatial_dimension(););
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityId identifier([[maybe_unused]] stk::mesh::Entity entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return bulk->identifier(entity);)
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityRank entity_rank([[maybe_unused]] stk::mesh::Entity entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return bulk->entity_rank(entity););
  }

  KOKKOS_FUNCTION
  stk::mesh::EntityKey entity_key([[maybe_unused]] stk::mesh::Entity entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return bulk->entity_key(entity););
  }

  KOKKOS_FUNCTION
  unsigned local_id([[maybe_unused]] stk::mesh::Entity entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return bulk->local_id(entity););
  }

  KOKKOS_FUNCTION
  stk::mesh::Entity get_entity([[maybe_unused]] stk::mesh::EntityRank rank,
                               [[maybe_unused]] const stk::mesh::FastMeshIndex& meshIndex) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return (*(bulk->buckets(rank)[meshIndex.bucket_id]))[meshIndex.bucket_ord];);
  }

  KOKKOS_FUNCTION
  stk::mesh::Entity linear_get_entity([[maybe_unused]] stk::mesh::EntityRank rank,
                                      [[maybe_unused]] stk::mesh::EntityId entityId) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return bulk->get_entity(rank, entityId););
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_connected_entities([[maybe_unused]] stk::mesh::EntityRank rank,
                                           [[maybe_unused]] const stk::mesh::FastMeshIndex &entity,
                                           [[maybe_unused]] stk::mesh::EntityRank connectedRank) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST((
      const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
      return bucket.get_connected_entities(entity.bucket_ord, connectedRank);
    ));
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_connected_ordinals([[maybe_unused]] stk::mesh::EntityRank rank,
                                           [[maybe_unused]] const stk::mesh::FastMeshIndex &entity,
                                           [[maybe_unused]] stk::mesh::EntityRank connectedRank) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST((
      const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
      return ConnectedOrdinals(bucket.begin_ordinals(entity.bucket_ord, connectedRank), bucket.num_connectivity(entity.bucket_ord, connectedRank));
    ));
  }

  KOKKOS_FUNCTION
  ConnectedNodes get_nodes([[maybe_unused]] stk::mesh::EntityRank rank,
                           [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_connected_entities(rank, entity, stk::topology::NODE_RANK););
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_edges([[maybe_unused]] stk::mesh::EntityRank rank,
                              [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_connected_entities(rank, entity, stk::topology::EDGE_RANK););
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_faces([[maybe_unused]] stk::mesh::EntityRank rank,
                              [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_connected_entities(rank, entity, stk::topology::FACE_RANK););
  }

  KOKKOS_FUNCTION
  ConnectedEntities get_elements([[maybe_unused]] stk::mesh::EntityRank rank,
                                 [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_connected_entities(rank, entity, stk::topology::ELEM_RANK););
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_node_ordinals([[maybe_unused]] stk::mesh::EntityRank rank,
                                      [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_connected_ordinals(rank, entity, stk::topology::NODE_RANK););
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_edge_ordinals([[maybe_unused]] stk::mesh::EntityRank rank,
                                      [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_connected_ordinals(rank, entity, stk::topology::EDGE_RANK););
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_face_ordinals([[maybe_unused]] stk::mesh::EntityRank rank,
                                      [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_connected_ordinals(rank, entity, stk::topology::FACE_RANK););
  }

  KOKKOS_FUNCTION
  ConnectedOrdinals get_element_ordinals([[maybe_unused]] stk::mesh::EntityRank rank,
                                         [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_connected_ordinals(rank, entity, stk::topology::ELEM_RANK););
  }

  KOKKOS_FUNCTION
  Permutations get_permutations([[maybe_unused]] stk::mesh::EntityRank rank,
                                [[maybe_unused]] const stk::mesh::FastMeshIndex &entity,
                                [[maybe_unused]] stk::mesh::EntityRank connectedRank) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST((
      const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
      return Permutations(bucket.begin_permutations(entity.bucket_ord, connectedRank), bucket.num_connectivity(entity.bucket_ord, connectedRank));
    ));
  }

  KOKKOS_FUNCTION
  Permutations get_node_permutations([[maybe_unused]] stk::mesh::EntityRank rank,
                                     [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_permutations(rank, entity, stk::topology::NODE_RANK););
  }

  KOKKOS_FUNCTION
  Permutations get_edge_permutations([[maybe_unused]] stk::mesh::EntityRank rank,
                                     [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_permutations(rank, entity, stk::topology::EDGE_RANK););
  }

  KOKKOS_FUNCTION
  Permutations get_face_permutations([[maybe_unused]] stk::mesh::EntityRank rank,
                                     [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_permutations(rank, entity, stk::topology::FACE_RANK););
  }

  KOKKOS_FUNCTION
  Permutations get_element_permutations([[maybe_unused]] stk::mesh::EntityRank rank,
                                        [[maybe_unused]] const stk::mesh::FastMeshIndex &entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return get_permutations(rank, entity, stk::topology::ELEM_RANK););
  }

  KOKKOS_FUNCTION
  stk::mesh::FastMeshIndex fast_mesh_index([[maybe_unused]] stk::mesh::Entity entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST((
      const stk::mesh::MeshIndex &meshIndex = bulk->mesh_index(entity);
      return stk::mesh::FastMeshIndex{meshIndex.bucket->bucket_id(), static_cast<unsigned>(meshIndex.bucket_ordinal)};
    ));
  }

  KOKKOS_FUNCTION
  stk::mesh::FastMeshIndex device_mesh_index([[maybe_unused]] stk::mesh::Entity entity) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return fast_mesh_index(entity););
  }

  KOKKOS_FUNCTION
  stk::NgpVector<unsigned> get_bucket_ids([[maybe_unused]] stk::mesh::EntityRank rank,
                                          [[maybe_unused]] const stk::mesh::Selector &selector) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return stk::mesh::get_bucket_ids(*bulk, rank, selector););
  }

  KOKKOS_FUNCTION
  unsigned num_buckets([[maybe_unused]] stk::mesh::EntityRank rank) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
    KOKKOS_IF_ON_HOST(return bulk->buckets(rank).size(););
  }

  KOKKOS_FUNCTION
  const BucketType & get_bucket([[maybe_unused]] stk::mesh::EntityRank rank,
                                [[maybe_unused]] unsigned i) const
  {
    KOKKOS_IF_ON_DEVICE((STK_NGP_ThrowErrorMsg("HostMesh only works on CPU/HOST.")));
#if !defined(STK_ENABLE_GPU) && !defined(NDEBUG)
    stk::mesh::EntityRank numRanks = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
    STK_NGP_ThrowAssert(rank < numRanks);
    STK_NGP_ThrowAssert(i < bulk->buckets(rank).size());
#endif
    KOKKOS_IF_ON_HOST((
    return *bulk->buckets(rank)[i];
    ));
  }

  NgpCommMapIndicesHostMirror<stk::ngp::MemSpace> volatile_fast_shared_comm_map(stk::topology::rank_t rank, int proc,
                                                                         bool includeGhosts=false) const
  {
    if (rank != cachedRank || proc != cachedProc) {
      cachedCommMap = bulk->template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(rank, proc, includeGhosts);
      cachedRank = rank;
      cachedProc = proc;
    }

    return cachedCommMap;
  }

  template <typename... EntityIdsParams, typename... AddPartParams, typename... EntitiesParams>
  void batch_declare_entities(stk::topology::rank_t rank,
                              const Kokkos::View<unsigned*, EntityIdsParams...>& entityIds,
                              const Kokkos::View<PartOrdinal*, AddPartParams...>& addPartOrdinals,
                              Kokkos::View<Entity*, EntitiesParams...>& requestedEntities)
  {
    using EntitiesMemorySpace = typename std::remove_reference<decltype(entityIds)>::type::memory_space;
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, EntitiesMemorySpace>::accessible,
                  "The memory space of the 'entities' View is inaccessible from the HostMesh execution space");

    Kokkos::resize(requestedEntities, entityIds.extent(0));
    std::vector<Entity> requestedEntitiesVector;

    std::vector<unsigned long> newIds(entityIds.extent(0));
    for (unsigned i = 0; i < newIds.size(); ++i) {
      newIds[i] = entityIds(i);
    }

    PartVector parts(addPartOrdinals.extent(0));
    for (unsigned i = 0; i < parts.size(); ++i) {
      parts[i] = bulk->mesh_meta_data().get_parts()[addPartOrdinals(i)];
    }

    bulk->modification_begin();
    bulk->declare_entities(rank, newIds, parts, requestedEntitiesVector);
    bulk->modification_end();
    for (unsigned i = 0; i < requestedEntities.extent(0); ++i) {
      requestedEntities(i) = requestedEntitiesVector[i];
    }
  }

  template <typename... EntitiesParams>
  void batch_destroy_entities(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities)
  {
    using EntitiesMemorySpace = typename std::remove_reference<decltype(entities)>::type::memory_space;
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, EntitiesMemorySpace>::accessible,
                  "The memory space of the 'entities' View is inaccessible from the HostMesh execution space");

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
                  "The memory space of the 'entities' View is inaccessible from the HostMesh execution space");
    using ResultsMemorySpace = typename std::remove_reference<decltype(wasDestroyed)>::type::memory_space;
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, ResultsMemorySpace>::accessible,
                  "The memory space of the 'wasDestroyed' View is inaccessible from the HostMesh execution space");
    STK_ThrowRequireMsg(wasDestroyed.extent(0) == entities.extent(0),
                        "batch_destroy_entities: 'wasDestroyed' View must be the same length as 'entities'");

    const unsigned numEntities = entities.extent(0);
    for (unsigned i = 0; i < numEntities; ++i) {
      wasDestroyed(i) = impl::can_destroy_entity(*bulk, entities(i));
    }

    bulk->modification_begin();
    for (unsigned i = 0; i < numEntities; ++i) {
      if (wasDestroyed(i)) {
        bulk->destroy_entity(entities(i));
      }
    }
    bulk->modification_end();
  }

  template <typename... EntitiesParams>
  void batch_destroy_relations(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                               stk::mesh::EntityRank connectedRank)
  {
    bulk->modification_begin();
    for (size_t i = 0; i < entities.extent(0); ++i) {
      stk::mesh::destroy_relations(*bulk, entities(i), connectedRank);
    }
    bulk->modification_end();
  }

  template <typename... FromParams, typename... OffsetParams, typename... ToParams,
            typename... OrdinalParams, typename... PermParams>
  void batch_declare_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
                               const Kokkos::View<unsigned*, OffsetParams...>& offsets,
                               const Kokkos::View<stk::mesh::Entity*, ToParams...>& toEntities,
                               const Kokkos::View<stk::mesh::RelationIdentifier*, OrdinalParams...>& ordinals,
                               const Kokkos::View<stk::mesh::Permutation*, PermParams...>& permutations)
  {
    impl_batch_declare_relations(fromEntities, offsets, toEntities, ordinals, permutations);
  }

  template <typename... FromParams, typename... OffsetParams, typename... ToParams,
            typename... OrdinalParams>
  void batch_declare_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
                               const Kokkos::View<unsigned*, OffsetParams...>& offsets,
                               const Kokkos::View<stk::mesh::Entity*, ToParams...>& toEntities,
                               const Kokkos::View<stk::mesh::RelationIdentifier*, OrdinalParams...>& ordinals)
  {
    Kokkos::View<stk::mesh::Permutation*, NgpMemSpace> permutations("invalid_permutations", toEntities.extent(0));
    for (size_t i = 0; i < permutations.extent(0); ++i) {
      permutations(i) = stk::mesh::Permutation::INVALID_PERMUTATION;
    }

    impl_batch_declare_relations(fromEntities, offsets, toEntities, ordinals, permutations);
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

  stk::mesh::BulkData &get_bulk_on_host()
  {
    return *bulk;
  }

  const stk::mesh::BulkData &get_bulk_on_host() const
  {
    return *bulk;
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
                  "The memory space of the 'entities' View is inaccessible from the HostMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, AddPartOrdinalsMemorySpace>::accessible,
                  "The memory space of the 'addPartOrdinals' View is inaccessible from the HostMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, RemovePartOrdinalsMemorySpace>::accessible,
                  "The memory space of the 'removePartOrdinals' View is inaccessible from the HostMesh execution space");

    std::vector<stk::mesh::Entity> hostEntities;
    std::vector<stk::mesh::Part*> hostAddParts;
    std::vector<stk::mesh::Part*> hostRemoveParts;

    hostEntities.reserve(entities.extent(0));
    for (size_t i = 0; i < entities.extent(0); ++i) {
      hostEntities.push_back(entities[i]);
    }

    const stk::mesh::PartVector& parts = bulk->mesh_meta_data().get_parts();

    hostAddParts.reserve(addPartOrdinals.extent(0));
    for (size_t i = 0; i < addPartOrdinals.extent(0); ++i) {
      const size_t partOrdinal = addPartOrdinals[i];
      STK_ThrowRequire(partOrdinal < parts.size());
      hostAddParts.push_back(parts[partOrdinal]);
    }

    hostRemoveParts.reserve(removePartOrdinals.extent(0));
    for (size_t i = 0; i < removePartOrdinals.extent(0); ++i) {
      const size_t partOrdinal = removePartOrdinals[i];
      STK_ThrowRequire(partOrdinal < parts.size());
      hostRemoveParts.push_back(parts[partOrdinal]);
    }

    bulk->batch_change_entity_parts(hostEntities, hostAddParts, hostRemoveParts);
  }

#ifndef STK_HIDE_DEPRECATED_CODE
  STK_DEPRECATED_MSG("Use update_bulk_data() instead.")
  void sync_to_host() {}

  STK_DEPRECATED_MSG("Use need_update_bulk_data() instead.")
  bool need_sync_to_host() const override {
    return false;
  }
#endif

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void impl_batch_change_entity_parts(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
  {
    batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  }

  auto& get_ngp_parallel_sum_host_buffer_offsets() {
    return impl::get_ngp_mesh_host_data<stk::ngp::MemSpace>(*bulk)->m_hostBufferOffsets;
  }

  auto& get_ngp_parallel_sum_host_mesh_indices_offsets() {
    return impl::get_ngp_mesh_host_data<stk::ngp::MemSpace>(*bulk)->m_hostMeshIndicesOffsets;
  }

  auto& get_ngp_parallel_sum_device_mesh_indices_offsets() {
    return impl::get_ngp_mesh_host_data<stk::ngp::MemSpace>(*bulk)->m_hostMeshIndicesOffsets;
  }

  auto& get_ngp_parallel_sum_host_byte_buffer() {
    return impl::get_ngp_mesh_host_data<stk::ngp::MemSpace>(*bulk)->m_byteBuffer;
  }

  auto& get_ngp_parallel_sum_device_byte_buffer() {
    return impl::get_ngp_mesh_host_data<stk::ngp::MemSpace>(*bulk)->m_byteBuffer;
  }

private:
  template <typename... FromParams, typename... OffsetParams, typename... ToParams,
            typename... OrdinalParams, typename... PermParams>
  void impl_batch_declare_relations(const Kokkos::View<stk::mesh::Entity*, FromParams...>& fromEntities,
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
                  "The memory space of the 'fromEntities' View is inaccessible from the HostMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, OffsetMemorySpace>::accessible,
                  "The memory space of the 'offsets' View is inaccessible from the HostMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, ToMemorySpace>::accessible,
                  "The memory space of the 'toEntities' View is inaccessible from the HostMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, OrdinalMemorySpace>::accessible,
                  "The memory space of the 'ordinals' View is inaccessible from the HostMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, PermutationsMemorySpace>::accessible,
                  "The memory space of the 'permutations' View is inaccessible from the HostMesh execution space");

    const size_t numFrom = fromEntities.extent(0);

#ifndef NDEBUG
    // Validate the CRS layout.
    STK_ThrowRequire(permutations.extent(0) == toEntities.extent(0));
    STK_ThrowRequire(ordinals.extent(0) == toEntities.extent(0));
    if (numFrom > 0) {
      STK_ThrowRequire(offsets.extent(0) == numFrom + 1);
      STK_ThrowRequire(offsets(0) == 0u);
      for (size_t i = 0; i < numFrom; ++i) {
        STK_ThrowRequire(offsets(i) <= offsets(i + 1));
      }
      STK_ThrowRequire(offsets(numFrom) == toEntities.extent(0));
    }
#endif

    stk::mesh::OrdinalVector scratch1, scratch2, scratch3;
    bulk->modification_begin();
    for (size_t i = 0; i < numFrom; ++i) {
      for (unsigned j = offsets(i); j < offsets(i + 1); ++j) {
        bulk->declare_relation(fromEntities(i), toEntities(j), ordinals(j),
                               permutations(j), scratch1, scratch2, scratch3);
      }
    }
    bulk->modification_end();
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
    using FromMemorySpace = typename std::remove_reference<decltype(fromEntities)>::type::memory_space;
    using ToMemorySpace = typename std::remove_reference<decltype(toEntities)>::type::memory_space;
    using PermMemorySpace = typename std::remove_reference<decltype(permutations)>::type::memory_space;

    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, FromMemorySpace>::accessible,
                  "The memory space of the 'fromEntities' View is inaccessible from the HostMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, ToMemorySpace>::accessible,
                  "The memory space of the 'toEntities' View is inaccessible from the HostMesh execution space");
    static_assert(Kokkos::SpaceAccessibility<MeshExecSpace, PermMemorySpace>::accessible,
                  "The memory space of the 'permutations' View is inaccessible from the HostMesh execution space");

    const size_t numFrom = fromEntities.extent(0);
    const size_t numCols = toEntities.extent(1);
    const bool hasPermutations = permutations.extent(0) > 0;

#ifndef NDEBUG
    if (numFrom > 0) {
      STK_ThrowRequire(toEntities.extent(0) == numFrom);   // one row of connectivity per from-entity
      STK_ThrowRequire(bulk->is_valid(fromEntities(0)));
      const stk::topology topo = bulk->bucket(fromEntities(0)).topology();
      const unsigned expectedCols = topo.num_sub_topology(connectedRank);
      STK_ThrowRequire(expectedCols > 0);          // connectedRank must be a valid downward sub-rank
      STK_ThrowRequire(numCols == expectedCols);   // full complement required
      for (size_t i = 0; i < numFrom; ++i) {
        STK_ThrowRequire(bulk->is_valid(fromEntities(i)));
        STK_ThrowRequire(bulk->bucket(fromEntities(i)).topology() == topo);  // single uniform topology
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
    for (size_t i = 0; i < numFrom; ++i) {
      const size_t base = i * numCols;
      offsets(i) = static_cast<unsigned>(base);
      for (size_t j = 0; j < numCols; ++j) {
        flatToEntities(base + j) = toEntities(i, j);
        flatOrdinals(base + j) = static_cast<stk::mesh::RelationIdentifier>(j);
        if (hasPermutations) {
          flatPermutations(base + j) = permutations(i, j);
        }
      }
    }
    if (numFrom > 0) {
      offsets(numFrom) = static_cast<unsigned>(numRelations);
    }
  }

  stk::mesh::BulkData *bulk;
  size_t m_syncCountWhenUpdated;
  mutable stk::mesh::EntityRank cachedRank = stk::topology::INVALID_RANK;
  mutable int cachedProc = -1;
  mutable NgpCommMapIndicesHostMirror<stk::ngp::MemSpace> cachedCommMap;
};

using HostMesh = HostMeshT<typename stk::ngp::HostExecSpace::memory_space>;

}
}

#endif
