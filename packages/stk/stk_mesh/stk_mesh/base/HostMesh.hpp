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

  HostMeshT()
    : NgpMeshBase(),
      bulk(nullptr)
  {

  }

  HostMeshT(const stk::mesh::BulkData& b)
    : NgpMeshBase(),
      bulk(&const_cast<stk::mesh::BulkData&>(b)),
      m_syncCountWhenUpdated(bulk->synchronized_count())
  {
    require_ngp_mesh_rank_limit(bulk->mesh_meta_data());
  }

  virtual ~HostMeshT() override = default;

  HostMeshT(const HostMeshT &) = default;
  HostMeshT(HostMeshT &&) = default;
  HostMeshT& operator=(const HostMeshT&) = default;
  HostMeshT& operator=(HostMeshT&&) = default;

  void update_mesh() override
  {
    m_syncCountWhenUpdated = bulk->synchronized_count();
  }

  unsigned get_spatial_dimension() const
  {
    return bulk->mesh_meta_data().spatial_dimension();
  }

  stk::mesh::EntityId identifier(stk::mesh::Entity entity) const
  {
    return bulk->identifier(entity);
  }

  stk::mesh::EntityRank entity_rank(stk::mesh::Entity entity) const
  {
    return bulk->entity_rank(entity);
  }

  stk::mesh::EntityKey entity_key(stk::mesh::Entity entity) const
  {
    return bulk->entity_key(entity);
  }

  stk::mesh::Entity get_entity(stk::mesh::EntityRank rank,
                               const stk::mesh::FastMeshIndex& meshIndex) const
  {
    return (*(bulk->buckets(rank)[meshIndex.bucket_id]))[meshIndex.bucket_ord];
  }

  ConnectedEntities get_connected_entities(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
  {
    const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
    return bucket.get_connected_entities(entity.bucket_ord, connectedRank);
  }

  ConnectedOrdinals get_connected_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
  {
    const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
    return ConnectedOrdinals(bucket.begin_ordinals(entity.bucket_ord, connectedRank), bucket.num_connectivity(entity.bucket_ord, connectedRank));
  }

  ConnectedNodes get_nodes(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::NODE_RANK);
  }

  ConnectedEntities get_edges(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::EDGE_RANK);
  }

  ConnectedEntities get_faces(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::FACE_RANK);
  }

  ConnectedEntities get_elements(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_entities(rank, entity, stk::topology::ELEM_RANK);
  }

  ConnectedOrdinals get_node_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::NODE_RANK);
  }

  ConnectedOrdinals get_edge_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::EDGE_RANK);
  }

  ConnectedOrdinals get_face_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::FACE_RANK);
  }

  ConnectedOrdinals get_element_ordinals(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_connected_ordinals(rank, entity, stk::topology::ELEM_RANK);
  }

  Permutations get_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity, stk::mesh::EntityRank connectedRank) const
  {
    const stk::mesh::Bucket& bucket = get_bucket(rank, entity.bucket_id);
    return Permutations(bucket.begin_permutations(entity.bucket_ord, connectedRank), bucket.num_connectivity(entity.bucket_ord, connectedRank));
  }

  Permutations get_node_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::NODE_RANK);
  }

  Permutations get_edge_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::EDGE_RANK);
  }

  Permutations get_face_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::FACE_RANK);
  }

  Permutations get_element_permutations(stk::mesh::EntityRank rank, const stk::mesh::FastMeshIndex &entity) const
  {
    return get_permutations(rank, entity, stk::topology::ELEM_RANK);
  }

  stk::mesh::FastMeshIndex fast_mesh_index(stk::mesh::Entity entity) const
  {
    const stk::mesh::MeshIndex &meshIndex = bulk->mesh_index(entity);
    return stk::mesh::FastMeshIndex{meshIndex.bucket->bucket_id(), static_cast<unsigned>(meshIndex.bucket_ordinal)};
  }

  stk::mesh::FastMeshIndex device_mesh_index(stk::mesh::Entity entity) const
  {
    return fast_mesh_index(entity);
  }

  stk::NgpVector<unsigned> get_bucket_ids(stk::mesh::EntityRank rank, const stk::mesh::Selector &selector) const
  {
    return stk::mesh::get_bucket_ids(*bulk, rank, selector);
  }

  unsigned num_buckets(stk::mesh::EntityRank rank) const
  {
    return bulk->buckets(rank).size();
  }

  const BucketType & get_bucket(stk::mesh::EntityRank rank, unsigned i) const
  {
#ifndef NDEBUG
    stk::mesh::EntityRank numRanks = static_cast<stk::mesh::EntityRank>(bulk->mesh_meta_data().entity_rank_count());
    STK_NGP_ThrowAssert(rank < numRanks);
    STK_NGP_ThrowAssert(i < bulk->buckets(rank).size());
#endif
    return *bulk->buckets(rank)[i];
  }

  NgpCommMapIndicesHostMirrorT<NgpMemSpace> volatile_fast_shared_comm_map(stk::topology::rank_t rank, int proc) const
  {
    return bulk->volatile_fast_shared_comm_map<NgpMemSpace>(rank, proc);
  }

  stk::mesh::BulkData &get_bulk_on_host()
  {
    return *bulk;
  }

  const stk::mesh::BulkData &get_bulk_on_host() const
  {
    return *bulk;
  }

  bool is_up_to_date() const
  {
    return m_syncCountWhenUpdated == bulk->synchronized_count();
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

  void sync_to_host() {}

  bool need_sync_to_host() const override {
    return false;
  }

  template <typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
  void impl_batch_change_entity_parts(const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                 const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
  {
    batch_change_entity_parts(entities, addPartOrdinals, removePartOrdinals);
  }

private:
  stk::mesh::BulkData *bulk;
  size_t m_syncCountWhenUpdated;
};

using HostMesh = HostMeshT<typename stk::ngp::HostExecSpace::memory_space>;

}
}

#endif
