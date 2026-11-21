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

#ifndef STK_MESH_NGPMESHENTITYLIFECYCLE_HPP
#define STK_MESH_NGPMESHENTITYLIFECYCLE_HPP

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/EntityKey.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/NgpViewAlgorithms.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "Kokkos_Core.hpp"
#include <cstddef>
#include <vector>

namespace stk::mesh::impl {

template <typename EntityKeyViewType, typename LocalIdViewType, typename ExecSpace, typename... EntitiesParams>
void invalidate_destroyed_entity_maps(EntityKeyViewType entityKeys, LocalIdViewType entityLocalIds,
                                      const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& destroyedEntities,
                                      ExecSpace const& execSpace)
{
  const unsigned invalidLocalId = stk::mesh::GetInvalidLocalId();
  Kokkos::parallel_for("invalidate_destroyed_entity_keys",
    Kokkos::RangePolicy<ExecSpace>(execSpace, 0, destroyedEntities.extent(0)),
    KOKKOS_LAMBDA(const int i) {
      const stk::mesh::Entity entity = destroyedEntities(i);
      const unsigned offset = entity.local_offset();
      entityKeys(offset) = stk::mesh::EntityKey();
      entityLocalIds(offset) = invalidLocalId;
    }
  );
  Kokkos::fence();
}

template <typename MeshType, typename ResultView, typename... EntitiesParams>
Kokkos::View<stk::mesh::Entity*, typename MeshType::ngp_mem_space>
get_destroyable_entities(const MeshType& ngpMesh,
                         const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                         const ResultView& wasDestroyed)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using NgpMemSpace = typename MeshType::ngp_mem_space;

  Kokkos::View<stk::mesh::Entity*, NgpMemSpace> candidates(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "destroyableEntities"), entities.extent(0));
  Kokkos::deep_copy(candidates, entities);

  const stk::mesh::EntityRank meshEndRank = ngpMesh.get_end_rank();
  Kokkos::parallel_for("mark_destroyable_entities",
    Kokkos::RangePolicy<ExecSpace>(0, candidates.extent(0)),
    KOKKOS_LAMBDA(const int i) {
      const stk::mesh::Entity entity = candidates(i);
      bool destroyable = entity.is_local_offset_valid();
      if (destroyable) {
        const stk::mesh::EntityRank rank = ngpMesh.entity_rank(entity);
        const stk::mesh::FastMeshIndex idx = ngpMesh.device_mesh_index(entity);
        for (stk::mesh::EntityRank conRank = static_cast<stk::mesh::EntityRank>(rank + 1);
             conRank < meshEndRank; ++conRank) {
          if (ngpMesh.get_connected_entities(rank, idx, conRank).size() > 0) {
            destroyable = false;
            break;
          }
        }
      }
      wasDestroyed(i) = destroyable;
      if (!destroyable) {
        candidates(i) = stk::mesh::Entity();
      }
    }
  );
  Kokkos::fence();

  remove_invalid_entities_sort_unique_and_resize(candidates, ExecSpace{});
  return candidates;
}

template <typename MeshType, typename... EntitiesParams>
stk::mesh::EntityRank get_max_entity_rank(const MeshType& ngpMesh,
                                          const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  int maxRank = stk::topology::BEGIN_RANK;
  if (entities.extent(0) > 0) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, entities.extent(0)),
      KOKKOS_LAMBDA(const int i, int& localMax) {
        const int r = static_cast<int>(ngpMesh.entity_rank(entities(i)));
        localMax = (r > localMax) ? r : localMax;
      }, Kokkos::Max<int>(maxRank)
    );
  }
  return static_cast<stk::mesh::EntityRank>(maxRank);
}

template <typename RequestedEntitiesView, typename EntityKeyViewType, typename LocalIdViewType, typename EntityIdsView>
void fill_requested_entities_and_keys(RequestedEntitiesView requestedEntities,
                                      EntityKeyViewType entityKeys,
                                      LocalIdViewType entityLocalIds,
                                      const EntityIdsView& entityIds,
                                      stk::topology::rank_t rank,
                                      unsigned numInitialEntities,
                                      unsigned entityKeysInitialIndex)
{
  Kokkos::parallel_for("Create requested entities and keys", entityIds.extent(0), KOKKOS_LAMBDA(size_t i) {
    requestedEntities(i) = stk::mesh::Entity(numInitialEntities + i + 1);
    entityKeys(entityKeysInitialIndex + i) = stk::mesh::EntityKey(rank, entityIds(i));
    auto startId = (numInitialEntities == 0) ? 0 : entityLocalIds(numInitialEntities - 1);
    entityLocalIds(numInitialEntities + i) = startId + i;
  });
}

template <typename MeshIndicesView, typename RequestedEntitiesView>
void invalidate_new_entity_mesh_indices(MeshIndicesView deviceMeshIndices,
                                        RequestedEntitiesView requestedEntities,
                                        size_t numNewEntities)
{
  Kokkos::parallel_for("Invalidate new entity mesh indices", numNewEntities,
    KOKKOS_LAMBDA(size_t i) {
      deviceMeshIndices(requestedEntities(i).local_offset()) =
          stk::mesh::FastMeshIndex{stk::mesh::INVALID_BUCKET_ID, stk::mesh::INVALID_BUCKET_ID};
    });
}

inline stk::topology get_topology_for_entities_with_nodes(const stk::mesh::MetaData& meta,
                                                          stk::topology::rank_t rank,
                                                          const std::vector<stk::mesh::PartOrdinal>& addPartOrdinals,
                                                          size_t numNodesPerEntity)
{
  STK_ThrowRequireMsg(rank > stk::topology::NODE_RANK && rank < meta.entity_rank_count(),
                      "batch_declare_entities with connected nodes: rank " << rank << " must be above NODE_RANK");

  stk::topology topo = stk::topology::INVALID_TOPOLOGY;
  for (stk::mesh::PartOrdinal ord : addPartOrdinals) {
    const stk::topology partTopo = meta.get_part(ord).topology();
    if (partTopo == stk::topology::INVALID_TOPOLOGY) { continue; }
    STK_ThrowRequireMsg(topo == stk::topology::INVALID_TOPOLOGY || topo == partTopo,
                        "batch_declare_entities with connected nodes: add parts carry conflicting topologies "
                        << topo << " and " << partTopo);
    topo = partTopo;
  }

  STK_ThrowRequireMsg(topo != stk::topology::INVALID_TOPOLOGY,
                      "batch_declare_entities with connected nodes: add parts must include a part with a topology");
  STK_ThrowRequireMsg(topo.rank() == rank,
                      "batch_declare_entities with connected nodes: topology " << topo << " has rank " << topo.rank()
                      << ", not the declared rank " << rank);
  STK_ThrowRequireMsg(topo.num_nodes() == numNodesPerEntity,
                      "batch_declare_entities with connected nodes: topology " << topo << " has " << topo.num_nodes()
                      << " nodes, but nodeIds has " << numNodesPerEntity << " columns");
  return topo;
}

} // namespace stk::mesh::impl

#endif
