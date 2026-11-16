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

#ifndef STK_MESH_NGPMESHCONNECTIVITY_HPP
#define STK_MESH_NGPMESHCONNECTIVITY_HPP

#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/NgpViewAlgorithms.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include <cstddef>

namespace stk::mesh::impl {

struct RelationToDestroy {
  stk::mesh::Entity from;
  stk::mesh::Entity to;
  stk::mesh::ConnectivityOrdinal ord = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;

  KOKKOS_INLINE_FUNCTION
  bool operator==(RelationToDestroy const& rhs) const {
    return from == rhs.from && to == rhs.to && ord == rhs.ord;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(RelationToDestroy const& rhs) const {
    return !(*this == rhs);
  }
};

struct RelationCompareByTo
{
  KOKKOS_INLINE_FUNCTION
  bool operator()(RelationToDestroy const& lhs, RelationToDestroy const& rhs) const {
    return lhs.to.local_offset() < rhs.to.local_offset();
  }
};

// The inverse (declare) operation carries the same {from, to, ord} triple as RelationToDestroy and
// sorts by to-entity the same way, so it reuses both directly rather than duplicating the type.
using RelationToDeclare = RelationToDestroy;

template <typename ExecSpace, typename EntityViewType>
void throw_on_invalid_input_entities(const EntityViewType& entities)
{
  unsigned numInvalidEntities = 0;
  Kokkos::parallel_reduce("check_valid_input_entities",
    Kokkos::RangePolicy<ExecSpace>(0, entities.extent(0)),
    KOKKOS_LAMBDA(const int i, unsigned& localCount) {
      if (!entities(i).is_local_offset_valid()) { ++localCount; }
    }, numInvalidEntities);
  Kokkos::fence();
  STK_ThrowRequireMsg(numInvalidEntities == 0,
    "batch mesh-modification called with " << numInvalidEntities
    << " invalid entit" << (numInvalidEntities == 1 ? "y" : "ies"));
}

template <typename MeshType, typename EntityViewType>
unsigned get_max_num_connectivity(const MeshType& ngpMesh, const EntityViewType& entities,
                                  EntityRank connectedRank)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  unsigned max = 0;
  Kokkos::parallel_reduce("max_num_connectivity", Kokkos::RangePolicy<ExecSpace>(0, entities.extent(0)),
    KOKKOS_LAMBDA(const int i, unsigned& lmax) {
      const Entity e = entities(i);
      if (e.is_local_offset_valid()) {
        const EntityRank entityRank = ngpMesh.entity_rank(e);
        if (entityRank != connectedRank) {
          const FastMeshIndex idx = ngpMesh.device_mesh_index(e);
          const unsigned num = ngpMesh.get_connected_entities(entityRank, idx, connectedRank).size();
          lmax = (num > lmax) ? num : lmax;
        }
      }
    }, Kokkos::Max<unsigned>(max)
  );
  Kokkos::fence();

  return max;
}

template <typename MeshType, typename EntityViewType, typename RelationViewType>
void fill_relations_to_destroy(const MeshType& ngpMesh, const EntityViewType& entities,
                               EntityRank connectedRank, unsigned maxNumConn,
                               RelationViewType& relations)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using RelationType = typename RelationViewType::value_type;

  Kokkos::parallel_for("fill_relations_to_destroy", Kokkos::RangePolicy<ExecSpace>(0, entities.extent(0)),
    KOKKOS_LAMBDA(const int i) {
      const Entity entity = entities(i);
      if (!entity.is_local_offset_valid()) { return; }

      const EntityRank entityRank = ngpMesh.entity_rank(entity);
      if (entityRank == connectedRank) { return; }
      const bool downward = entityRank > connectedRank;

      const FastMeshIndex idx = ngpMesh.device_mesh_index(entity);
      auto connEntities = ngpMesh.get_connected_entities(entityRank, idx, connectedRank);
      auto connOrdinals = ngpMesh.get_connected_ordinals(entityRank, idx, connectedRank);

      const unsigned base = i * maxNumConn;
      for (unsigned connIdx = 0; connIdx < connEntities.size(); ++connIdx) {
        const Entity other = connEntities[connIdx];
        const ConnectivityOrdinal ord = connOrdinals[connIdx];
        relations(base + connIdx) = downward ? RelationType{entity, other, ord}
                                             : RelationType{other, entity, ord};
      }
    }
  );
  Kokkos::fence();
}

// Performance note: thread-per-from-entity strides the CRS reads (toEntities/ordinals) across a
// warp rather than coalescing them, but that is the right choice here.  This kernel is latency-
// bound on the per-from-entity connectivity gathers (which are uncoalesced under any policy), the
// CRS arrays are the minority of the traffic, monotonic offsets keep a warp's strided reads inside
// a small contiguous region that L1/L2 captures, and the slices are short (single-digit relations
// per entity).  A team-per-entity would coalesce the CRS reads but waste most of a warp's lanes on
// those short slices and shed occupancy -- the classic CSR scalar-vs-vector SpMV trade-off, where
// scalar (thread-per-row) wins for short rows.  If this ever profiles as hot and memory-bound, the
// coalescing-optimal (and load-balanced) evolution is NOT a team loop but a flat-over-relations
// RangePolicy(0, numRelations) -- one thread per flat j, so consecutive threads read consecutive
// toEntities(j)/ordinals(j).  That requires reformulating the intra-slice dedup to the mask-
// independent "is there an earlier equal (to, ord) entry in this slice" rule (dropping the
// sequential keep(p) read below), plus a j->from-entity back-map from the offsets.
template <typename MeshType, typename FromViewType, typename OffsetViewType,
          typename ToViewType, typename OrdinalViewType, typename MaskViewType>
void compute_new_relation_keep_mask(const MeshType& ngpMesh,
                                    const FromViewType& fromEntities,
                                    const OffsetViewType& offsets,
                                    const ToViewType& toEntities,
                                    const OrdinalViewType& ordinals,
                                    MaskViewType& keep)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  Kokkos::parallel_for("compute_new_relation_keep_mask",
    Kokkos::RangePolicy<ExecSpace>(0, fromEntities.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      const Entity from = fromEntities(i);
      const EntityRank fromRank = ngpMesh.entity_rank(from);
      const FastMeshIndex fromIdx = ngpMesh.device_mesh_index(from);

      const unsigned begin = offsets(i);
      const unsigned end = offsets(i + 1);
      for (unsigned j = begin; j < end; ++j) {
        const Entity to = toEntities(j);
        const EntityRank toRank = ngpMesh.entity_rank(to);
        const ConnectivityOrdinal ord = static_cast<ConnectivityOrdinal>(ordinals(j));

        bool duplicate = false;

        auto connEntities = ngpMesh.get_connected_entities(fromRank, fromIdx, toRank);
        auto connOrdinals = ngpMesh.get_connected_ordinals(fromRank, fromIdx, toRank);
        for (unsigned connIdx = 0; connIdx < connEntities.size(); ++connIdx) {
          if (connEntities[connIdx] == to && connOrdinals[connIdx] == ord) { duplicate = true; break; }
        }

        for (unsigned priorIdx = begin; !duplicate && priorIdx < j; ++priorIdx) {
          if (keep(priorIdx) && toEntities(priorIdx) == to &&
              static_cast<ConnectivityOrdinal>(ordinals(priorIdx)) == ord) {
            duplicate = true;
          }
        }

        keep(j) = !duplicate;
      }
    }
  );
  Kokkos::fence();
}

template <typename FromViewType, typename OffsetViewType, typename ToViewType,
          typename OrdinalViewType, typename MaskViewType, typename RelationViewType>
void fill_relations_to_declare(const FromViewType& fromEntities, const OffsetViewType& offsets,
                               const ToViewType& toEntities, const OrdinalViewType& ordinals,
                               const MaskViewType& keep, RelationViewType& relations)
{
  using ExecSpace = typename RelationViewType::execution_space;
  Kokkos::parallel_for("fill_relations_to_declare", Kokkos::RangePolicy<ExecSpace>(0, fromEntities.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      const Entity from = fromEntities(i);
      const unsigned begin = offsets(i);
      const unsigned end = offsets(i + 1);
      for (unsigned j = begin; j < end; ++j) {
        relations(j) = keep(j)
          ? RelationToDeclare{from, toEntities(j), static_cast<ConnectivityOrdinal>(ordinals(j))}
          : RelationToDeclare{};
      }
    }
  );
  Kokkos::fence();
}

template <typename RelationViewType, typename EntityViewType>
void collect_relation_to_entities(const RelationViewType& relations, EntityViewType& toEntities)
{
  using ExecSpace = typename EntityViewType::execution_space;
  Kokkos::parallel_for("collect_relation_to_entities", Kokkos::RangePolicy<ExecSpace>(0, relations.extent(0)),
    KOKKOS_LAMBDA(const size_t i) { toEntities(i) = relations(i).to; }
  );
  Kokkos::fence();
}

template <typename ExecSpace, typename EntityViewType>
unsigned compute_max_fan_in(EntityViewType const& toEntitiesWithDuplicates)
{
  const size_t n = toEntitiesWithDuplicates.extent(0);
  if (n == 0) { return 0; }

  STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, toEntitiesWithDuplicates));

  unsigned maxRun = 0;
  Kokkos::parallel_reduce("compute_max_fan_in", Kokkos::RangePolicy<ExecSpace>(0, n),
    KOKKOS_LAMBDA(const size_t i, unsigned& localMax) {
      if (i == 0 || !(toEntitiesWithDuplicates(i) == toEntitiesWithDuplicates(i - 1))) {
        unsigned run = 1;
        for (size_t k = i + 1; k < n && toEntitiesWithDuplicates(k) == toEntitiesWithDuplicates(i); ++k) {
          ++run;
        }
        if (run > localMax) { localMax = run; }
      }
    }, Kokkos::Max<unsigned>(maxRun));
  return maxRun;
}

template <typename MeshType, typename RelationViewType>
void validate_and_mark_relations(const MeshType& ngpMesh, RelationViewType& relations)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using RelationType = typename RelationViewType::value_type;

  Kokkos::parallel_for("validate_relations_to_destroy",
    Kokkos::RangePolicy<ExecSpace>(0, relations.extent(0)),
    KOKKOS_LAMBDA(const int i) {
      auto& rel = relations(i);
      const Entity fromEntity = rel.from;
      const Entity toEntity = rel.to;
      const ConnectivityOrdinal ord = rel.ord;

      bool valid = fromEntity.is_local_offset_valid() && toEntity.is_local_offset_valid();

      if (valid) {
        const EntityRank fromRank = ngpMesh.entity_rank(fromEntity);
        const EntityRank toRank = ngpMesh.entity_rank(toEntity);
        valid = (fromRank > toRank);

        if (valid) {
          const FastMeshIndex fromIdx = ngpMesh.device_mesh_index(fromEntity);
          auto connEntities = ngpMesh.get_connected_entities(fromRank, fromIdx, toRank);
          auto connOrdinals = ngpMesh.get_connected_ordinals(fromRank, fromIdx, toRank);

          bool found = false;
          for (unsigned connIdx = 0; connIdx < connEntities.size(); ++connIdx) {
            if (connOrdinals[connIdx] == ord && connEntities[connIdx] == toEntity) {
              found = true;
              break;
            }
          }
          valid = found;
        }
      }

      if (!valid) {
        rel = RelationType{};
      }
    }
  );
  Kokkos::fence();
}

struct DirectedConnectivityRemoval {
  stk::mesh::Entity owner;
  stk::mesh::Entity target;
  stk::mesh::EntityRank connRank = stk::mesh::InvalidEntityRank;
  stk::mesh::ConnectivityOrdinal ord = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;

  KOKKOS_INLINE_FUNCTION
  bool operator<(DirectedConnectivityRemoval const& rhs) const {
    return owner.local_offset() < rhs.owner.local_offset();
  }
};

// One directed connectivity addition: `owner` is the entity whose connectivity array grows,
// `target` the entity it connects to at `connRank`.  Each kept relation expands into two of
// these (downward from->to and the reciprocal upward to->from), exactly like the removal path.
// The sort key is composite -- (owner, connRank, ord, target) -- so that each (owner, connRank)
// run arrives in the same canonical order the host BucketConnDynamic keeps its slices in
// (ascending ordinal, then ascending entity; see find_sorted_insertion_index).  That lets
// insert_connectivities_no_grow merge a run straight into the existing slice and preserve the
// invariant.  Sorting on (ord, target) rather than an input sequence number also makes the key a
// total order on distinct entries, so Kokkos::sort not being stable cannot affect the result.
struct DirectedConnectivityAddition {
  stk::mesh::Entity owner;
  stk::mesh::Entity target;
  stk::mesh::EntityRank connRank = stk::mesh::InvalidEntityRank;
  stk::mesh::ConnectivityOrdinal ord = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
  stk::mesh::Permutation perm = stk::mesh::INVALID_PERMUTATION;

  KOKKOS_INLINE_FUNCTION
  bool operator<(DirectedConnectivityAddition const& rhs) const {
    if (owner.local_offset() != rhs.owner.local_offset()) {
      return owner.local_offset() < rhs.owner.local_offset();
    }
    if (connRank != rhs.connRank) {
      return connRank < rhs.connRank;
    }
    if (ord != rhs.ord) {
      return ord < rhs.ord;
    }
    return target.local_offset() < rhs.target.local_offset();
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(DirectedConnectivityAddition const& rhs) const {
    return owner == rhs.owner && target == rhs.target &&
           connRank == rhs.connRank && ord == rhs.ord;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(DirectedConnectivityAddition const& rhs) const {
    return !(*this == rhs);
  }
};

struct OwnerAdditionSizing {
  stk::mesh::Entity owner;
  unsigned numNew = 0;
  bool needPermutations = false;

  KOKKOS_INLINE_FUNCTION
  bool operator==(OwnerAdditionSizing const& rhs) const {
    return owner == rhs.owner && numNew == rhs.numNew && needPermutations == rhs.needPermutations;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(OwnerAdditionSizing const& rhs) const {
    return !(*this == rhs);
  }
};

template <typename MeshType, typename FromViewType, typename OffsetViewType, typename ToViewType,
          typename OrdinalViewType, typename PermViewType, typename MaskViewType, typename AdditionViewType>
void build_directed_additions(const MeshType& ngpMesh,
                              const FromViewType& fromEntities,
                              const OffsetViewType& offsets,
                              const ToViewType& toEntities,
                              const OrdinalViewType& ordinals,
                              const PermViewType& permutations,
                              bool hasPerm,
                              const MaskViewType& keep,
                              AdditionViewType& additions)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  Kokkos::parallel_for("build_directed_additions",
    Kokkos::RangePolicy<ExecSpace>(0, fromEntities.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      const Entity from = fromEntities(i);
      const EntityRank fromRank = ngpMesh.entity_rank(from);
      const unsigned begin = offsets(i);
      const unsigned end = offsets(i + 1);
      for (unsigned j = begin; j < end; ++j) {
        DirectedConnectivityAddition down{};
        DirectedConnectivityAddition up{};
        if (keep(j)) {
          const Entity to = toEntities(j);
          const EntityRank toRank = ngpMesh.entity_rank(to);
          const ConnectivityOrdinal ord = static_cast<ConnectivityOrdinal>(ordinals(j));
          const stk::mesh::Permutation perm =
              hasPerm ? permutations(j) : stk::mesh::INVALID_PERMUTATION;
          down = DirectedConnectivityAddition{from, to, toRank, ord, perm};
          up   = DirectedConnectivityAddition{to, from, fromRank, ord, perm};
        }
        additions(2 * j)     = down;
        additions(2 * j + 1) = up;
      }
    }
  );
  Kokkos::fence();
}

template <typename ExecSpace, typename AdditionViewType, typename SizingViewType>
void compute_owner_add_sizes(const AdditionViewType& sortedAdditions, SizingViewType& sizing)
{
  const unsigned numAdds = sortedAdditions.extent(0);
  Kokkos::parallel_for("compute_owner_add_sizes", Kokkos::RangePolicy<ExecSpace>(0, numAdds),
    KOKKOS_LAMBDA(const unsigned i) {
      if (i != 0 && sortedAdditions(i).owner == sortedAdditions(i - 1).owner) {
        sizing(i) = OwnerAdditionSizing{};
        return;
      }
      const stk::mesh::Entity owner = sortedAdditions(i).owner;
      unsigned count = 0;
      bool needPerm = false;
      unsigned j = i;
      while (j < numAdds && sortedAdditions(j).owner == owner) {
        needPerm = needPerm || stk::mesh::does_rank_have_valid_permutations(sortedAdditions(j).connRank);
        ++count;
        ++j;
      }
      sizing(i) = OwnerAdditionSizing{owner, count, needPerm};
    }
  );
  Kokkos::fence();
}

template <typename MeshType, typename AdditionViewType>
void insert_directed_additions_by_owner(const MeshType& ngpMesh, const AdditionViewType& sortedAdditions)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  const unsigned numAdds = sortedAdditions.extent(0);
  if (numAdds == 0) { return; }

  Kokkos::parallel_for("declare_by_owner", Kokkos::RangePolicy<ExecSpace>(0, numAdds),
    KOKKOS_LAMBDA(const unsigned i) {
      if (i != 0 && sortedAdditions(i).owner == sortedAdditions(i - 1).owner) { return; }
      auto& meshConn = ngpMesh.get_mesh_connectivity();
      const stk::mesh::Entity owner = sortedAdditions(i).owner;
      unsigned j = i;
      while (j < numAdds && sortedAdditions(j).owner == owner) { ++j; }
      meshConn.insert_connectivities_no_grow(owner, &sortedAdditions(i), j - i);
    }
  );
  Kokkos::fence();
}

template <typename MeshType, typename RelationViewType>
unsigned get_max_from_parts_for_induction(const MeshType& ngpMesh, const RelationViewType& relations)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using BucketType = typename MeshType::BucketType;

  unsigned maxParts = 0;
  Kokkos::parallel_reduce("get_max_from_parts_for_induction",
    Kokkos::RangePolicy<ExecSpace>(0, relations.extent(0)),
    KOKKOS_LAMBDA(const unsigned i, unsigned& lmax) {
      const auto rel = relations(i);
      const EntityRank fromRank = ngpMesh.entity_rank(rel.from);
      const EntityRank toRank = ngpMesh.entity_rank(rel.to);
      if (fromRank > toRank) {
        const FastMeshIndex fromIdx = ngpMesh.device_mesh_index(rel.from);
        const BucketType& bucket = ngpMesh.get_bucket(fromRank, fromIdx.bucket_id);
        auto partOrdinalsPair = bucket.superset_part_ordinals();
        const unsigned numParts = static_cast<unsigned>(partOrdinalsPair.second - partOrdinalsPair.first);
        lmax = (numParts > lmax) ? numParts : lmax;
      }
    }, Kokkos::Max<unsigned>(maxParts)
  );
  Kokkos::fence();

  return maxParts;
}

template <typename MeshType, typename RelationViewType>
void sever_relations(const MeshType& ngpMesh, const RelationViewType& relations)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using MemSpace = typename RelationViewType::memory_space;
  using RemovalViewType = Kokkos::View<DirectedConnectivityRemoval*, MemSpace>;

  const unsigned numRelations = relations.extent(0);
  if (numRelations == 0) { return; }

  // Expand each relation into its two directed removals.
  RemovalViewType removals(Kokkos::view_alloc(Kokkos::WithoutInitializing, "directedRemovals"), numRelations * 2);
  Kokkos::parallel_for("build_directed_removals", Kokkos::RangePolicy<ExecSpace>(0, numRelations),
    KOKKOS_LAMBDA(const unsigned i) {
      const auto rel = relations(i);
      const EntityRank fromRank = ngpMesh.entity_rank(rel.from);
      const EntityRank toRank = ngpMesh.entity_rank(rel.to);
      removals(2 * i)     = DirectedConnectivityRemoval{rel.from, rel.to,   toRank,   rel.ord};
      removals(2 * i + 1) = DirectedConnectivityRemoval{rel.to,   rel.from, fromRank, rel.ord};
    }
  );
  Kokkos::fence();

  // Group removals so each owner's array is contiguous.
  Kokkos::sort(removals);

  const unsigned numRemovals = removals.extent(0);
  Kokkos::parallel_for("sever_by_owner", Kokkos::RangePolicy<ExecSpace>(0, numRemovals),
    KOKKOS_LAMBDA(const unsigned i) {
      if (i != 0 && removals(i).owner == removals(i - 1).owner) { return; }
      auto& meshConn = ngpMesh.get_mesh_connectivity();
      const stk::mesh::Entity owner = removals(i).owner;
      // This owner's removals are contiguous [i, j); remove them all in one compaction pass.
      unsigned j = i;
      while (j < numRemovals && removals(j).owner == owner) { ++j; }
      meshConn.remove_connectivities(owner, &removals(i), j - i);
    }
  );
  Kokkos::fence();
}

template<typename RelationValueType, typename MemSpace>
struct SortedRelationsWithOffsets {
  Kokkos::View<RelationValueType*, MemSpace> sortedRelations;
  Kokkos::View<unsigned*, MemSpace> relationOffsets;
};

template<typename ExecSpace, typename EntityViewType, typename RelationViewType>
SortedRelationsWithOffsets<typename RelationViewType::value_type, typename RelationViewType::memory_space>
build_relation_offsets_by_to_entity(const EntityViewType& entities, const RelationViewType& relations)
{
  using RelationValueType = typename RelationViewType::value_type;
  using MemSpace = typename RelationViewType::memory_space;

  const unsigned numEntities = entities.extent(0);
  const unsigned numRelations = relations.extent(0);

  Kokkos::View<RelationValueType*, MemSpace> sortedRelations(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "sortedRelationsByTo"), numRelations);
  Kokkos::deep_copy(sortedRelations, relations);
  Kokkos::sort(sortedRelations, RelationCompareByTo{});

  Kokkos::View<unsigned*, MemSpace> relationOffsets(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "relationOffsetsByEntity"), numEntities + 1);
  Kokkos::parallel_for("build_relation_offsets_by_entity", Kokkos::RangePolicy<ExecSpace>(0, numEntities),
    KOKKOS_LAMBDA(const unsigned k) {
      const auto target = entities(k).local_offset();
      unsigned lo = 0;
      unsigned hi = numRelations;
      while (lo < hi) {
        const unsigned mid = lo + (hi - lo) / 2;
        if (sortedRelations(mid).to.local_offset() < target) { lo = mid + 1; }
        else { hi = mid; }
      }
      relationOffsets(k) = lo;
      if (k == numEntities - 1) { relationOffsets(numEntities) = numRelations; }
    }
  );
  Kokkos::fence();

  return {sortedRelations, relationOffsets};
}

template<typename MeshType, typename EntityViewType, typename RelationViewType,
         typename NewPartsViewType, typename PartOrdinalsProxyViewType>
void set_new_part_list_per_entity_after_relation_removal(const MeshType& ngpMesh,
                                                         const EntityViewType& entities,
                                                         const RelationViewType& relations,
                                                         unsigned maxPartsPerEntity,
                                                         NewPartsViewType& newPartsPerEntity,
                                                         PartOrdinalsProxyViewType& partOrdinalsProxyView)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using TeamMember = typename stk::ngp::TeamPolicy<ExecSpace>::member_type;

  const unsigned numEntities = entities.extent(0);

  auto relOffsets = build_relation_offsets_by_to_entity<ExecSpace>(entities, relations);
  auto sortedRelations = relOffsets.sortedRelations;
  auto relationOffsets = relOffsets.relationOffsets;

  auto teamPolicy = stk::ngp::TeamPolicy<ExecSpace>(numEntities, Kokkos::AUTO);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for("set_new_part_lists_after_relation_removal", teamPolicy,
    KOKKOS_LAMBDA(const TeamMember& team) {
      auto idx = team.league_rank();
      auto myStartIdx = idx * maxPartsPerEntity;
      auto rank = ngpMesh.entity_rank(entities(idx));
      auto fastMeshIdx = ngpMesh.device_mesh_index(entities(idx));
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIdx.bucket_id);
      auto currentPartOrdsPair = bucket.superset_part_ordinals();

      Ordinal* dest = &newPartsPerEntity(myStartIdx);
      const Ordinal* first = currentPartOrdsPair.first;
      const Ordinal* last = currentPartOrdsPair.second;
      unsigned myNumParts = last - first;

      Kokkos::single(Kokkos::PerTeam(team),[&]() {
        my_copy(first, last, dest);

        typename PartOrdinalsProxyViewType::value_type proxyIndices(rank, dest, myNumParts);
        partOrdinalsProxyView(idx) = proxyIndices;
      });

      team.team_barrier();

      auto isRemovedPart = [&](Ordinal partOrdinal) { return partOrdinal == InvalidPartOrdinal; };

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, myNumParts),
        [&](const int& thIdx) {
          auto myPartOrdinalIdx = myStartIdx + thIdx;
          auto partToCheck = newPartsPerEntity(myPartOrdinalIdx);
          auto partToCheckRank = deviceBucketRepo.get_part_rank(partToCheck);
          auto inducingPart = deviceBucketRepo.is_inducible_part(partToCheck);

          if (!inducingPart || partToCheckRank <= rank) {
            return;
          }

          bool inducedBySeveredFrom = false;
          const unsigned relBegin = relationOffsets(idx);
          const unsigned relEnd = relationOffsets(idx + 1);
          for (unsigned r = relBegin; r < relEnd; ++r) {
            // All relations in [relBegin, relEnd) already have to == entities(idx).
            auto fromEntity = sortedRelations(r).from;
            if (ngpMesh.entity_rank(fromEntity) != partToCheckRank) { continue; }
            auto fromMeshIdx = ngpMesh.device_mesh_index(fromEntity);
            auto& fromBucket = ngpMesh.get_bucket(partToCheckRank, fromMeshIdx.bucket_id);
            if (fromBucket.member(partToCheck)) {
              inducedBySeveredFrom = true;
              break;
            }
          }

          if (!inducedBySeveredFrom) {
            return;
          }

          auto connectedUpperRankEntities = ngpMesh.get_connected_entities(rank, fastMeshIdx, partToCheckRank);
          bool stillInduced = false;

          Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, connectedUpperRankEntities.size()),
            [&](const int& connIdx, bool& localInduced) {
              auto entity = connectedUpperRankEntities[connIdx];
              auto bucketId = ngpMesh.device_mesh_index(entity).bucket_id;
              auto& localBucket = ngpMesh.get_bucket(partToCheckRank, bucketId);
              localInduced |= localBucket.member(partToCheck);
            }, Kokkos::LOr<bool>(stillInduced)
          );

          Kokkos::single(Kokkos::PerThread(team), [&]() {
            if (!stillInduced) {
              newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
            }
          });
        }
      );

      auto begin = Kokkos::Experimental::begin(newPartsPerEntity) + myStartIdx;
      auto end = begin + myNumParts;

      end = Kokkos::Experimental::remove_if(team, begin, end, isRemovedPart);
      end = Kokkos::Experimental::unique(team, begin, end);

      Kokkos::single(Kokkos::PerTeam(team),[&]() {
        partOrdinalsProxyView(idx).length = static_cast<unsigned>(Kokkos::Experimental::distance(begin, end));
      });
    }
  );

  Kokkos::fence();
}

template<typename MeshType, typename EntityViewType, typename RelationViewType,
         typename NewPartsViewType, typename PartOrdinalsProxyViewType>
void set_new_part_list_per_entity_after_relation_addition(const MeshType& ngpMesh,
                                                          const EntityViewType& entities,
                                                          const RelationViewType& relations,
                                                          unsigned maxPartsPerEntity,
                                                          NewPartsViewType& newPartsPerEntity,
                                                          PartOrdinalsProxyViewType& partOrdinalsProxyView)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  const unsigned numEntities = entities.extent(0);

  auto relOffsets = build_relation_offsets_by_to_entity<ExecSpace>(entities, relations);
  auto sortedRelations = relOffsets.sortedRelations;
  auto relationOffsets = relOffsets.relationOffsets;

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for("set_new_part_lists_after_relation_addition",
    Kokkos::RangePolicy<ExecSpace>(0, numEntities),
    KOKKOS_LAMBDA(const unsigned idx) {
      const auto myStartIdx = idx * maxPartsPerEntity;
      const auto rank = ngpMesh.entity_rank(entities(idx));
      const auto fastMeshIdx = ngpMesh.device_mesh_index(entities(idx));
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIdx.bucket_id);
      auto currentPartOrdsPair = bucket.superset_part_ordinals();

      Ordinal* dest = &newPartsPerEntity(myStartIdx);
      const Ordinal* first = currentPartOrdsPair.first;
      const Ordinal* last = currentPartOrdsPair.second;

      unsigned numParts = static_cast<unsigned>(my_copy(first, last, dest) - dest);

      const unsigned relBegin = relationOffsets(idx);
      const unsigned relEnd = relationOffsets(idx + 1);
      for (unsigned relIdx = relBegin; relIdx < relEnd; ++relIdx) {
        const auto fromEntity = sortedRelations(relIdx).from;
        const auto fromRank = ngpMesh.entity_rank(fromEntity);
        if (fromRank <= rank) { continue; }  // only higher-rank entities induce downward
        const auto fromMeshIdx = ngpMesh.device_mesh_index(fromEntity);
        auto& fromBucket = ngpMesh.get_bucket(fromRank, fromMeshIdx.bucket_id);
        auto fromPartOrdsPair = fromBucket.superset_part_ordinals();
        for (const Ordinal* fromPartOrd = fromPartOrdsPair.first; fromPartOrd != fromPartOrdsPair.second; ++fromPartOrd) {
          const Ordinal partOrdinal = *fromPartOrd;
          if (deviceBucketRepo.is_inducible_part(partOrdinal) &&
              deviceBucketRepo.get_part_rank(partOrdinal) == fromRank) {
            dest[numParts++] = partOrdinal;
          }
        }
      }

      // Insertion sort the collected part ordinals in place.
      for (unsigned sortIdx = 1; sortIdx < numParts; ++sortIdx) {
        const Ordinal key = dest[sortIdx];
        int priorIdx = static_cast<int>(sortIdx) - 1;
        while (priorIdx >= 0 && dest[priorIdx] > key) { dest[priorIdx + 1] = dest[priorIdx]; --priorIdx; }
        dest[priorIdx + 1] = key;
      }
      unsigned numUniqueParts = 0;
      for (unsigned partIdx = 0; partIdx < numParts; ++partIdx) {
        if (partIdx == 0 || dest[partIdx] != dest[numUniqueParts - 1]) { dest[numUniqueParts++] = dest[partIdx]; }
      }

      typename PartOrdinalsProxyViewType::value_type proxyIndices(rank, dest, numUniqueParts);
      partOrdinalsProxyView(idx) = proxyIndices;
    }
  );

  Kokkos::fence();
}

} // namespace stk::mesh::impl

#endif
