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

#ifndef STK_MESH_NGPMESHPARALLELCOMM_HPP
#define STK_MESH_NGPMESHPARALLELCOMM_HPP

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/EntityKey.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/NgpMeshPartChange.hpp"
#include "stk_mesh/baseImpl/NgpMeshTypes.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/stk_config.h"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/parallel/CommBuffer.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace stk::mesh::impl {

template <typename MeshType, typename CommMapSliceType>
KOKKOS_INLINE_FUNCTION
Entity find_entity_key_in_sorted_comm_map_range(MeshType const& ngpMesh, EntityRank rank,
                                                CommMapSliceType const& commMapSlice,
                                                size_t begin, size_t end, EntityKey key)
{
  size_t lo = begin;
  size_t hi = end;
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    auto midEntity = ngpMesh.get_entity(rank, commMapSlice(mid));
    if (ngpMesh.entity_key(midEntity) < key) { lo = mid + 1; } else { hi = mid; }
  }

  if (lo < end) {
    auto foundEntity = ngpMesh.get_entity(rank, commMapSlice(lo));
    if (ngpMesh.entity_key(foundEntity) == key) { return foundEntity; }
  }
  return Entity{};
}

template <typename MeshType>
KOKKOS_INLINE_FUNCTION
Entity find_entity_key_in_comm_map(MeshType const& ngpMesh, EntityRank rank, int proc, EntityKey key,
                                   bool includeGhosts)
{
  auto commMapSlice = ngpMesh.volatile_fast_shared_comm_map(static_cast<stk::topology::rank_t>(rank),
                                                            proc, includeGhosts);
  auto sharedSlice = ngpMesh.volatile_fast_shared_comm_map(static_cast<stk::topology::rank_t>(rank),
                                                           proc, false);
  auto numShared = sharedSlice.extent(0);

  auto foundEntity = find_entity_key_in_sorted_comm_map_range(ngpMesh, rank, commMapSlice, 0, numShared, key);
  if (foundEntity.is_local_offset_valid() || !includeGhosts) { return foundEntity; }

  return find_entity_key_in_sorted_comm_map_range(ngpMesh, rank, commMapSlice, numShared,
                                                  commMapSlice.extent(0), key);
}

template <typename EntityKeyWrapperViewType, typename ExecSpace>
void remove_invalid_entity_keys_sort_unique_and_resize(EntityKeyWrapperViewType& entityKeys, ExecSpace const& execSpace)
{
  using EntityKeyWrapper = typename EntityKeyWrapperViewType::value_type;
  using MemorySpace = typename EntityKeyWrapperViewType::memory_space;
  using EntityKeyWrapperUViewType = Kokkos::View<EntityKeyWrapper*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  auto isInvalidEntityKey = KOKKOS_LAMBDA(EntityKeyWrapper entityKeyWrapper) {
                              return !entityKeyWrapper.value.is_valid();
                            };
  auto end = Kokkos::Experimental::remove_if(execSpace, entityKeys, isInvalidEntityKey);

  auto begin = Kokkos::Experimental::begin(entityKeys);
  auto length = Kokkos::Experimental::distance(begin, end);

  EntityKeyWrapperUViewType uview(entityKeys.data(), length);
  Kokkos::sort(uview);
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, uview));

  merge_sorted_duplicate_wrapper_flags(uview, static_cast<size_t>(length), execSpace);

  auto newEnd = Kokkos::Experimental::unique(execSpace, uview);
  auto newBegin = Kokkos::Experimental::begin(uview);
  size_t newSize = Kokkos::Experimental::distance(newBegin, newEnd);
  if (newSize != entityKeys.extent(0)) {
    Kokkos::resize(entityKeys, newSize);
  }
}

template <typename MeshType, typename WrappedEntityViewType>
Kokkos::View<unsigned*> intersect_wrapped_entities_with_comm_map(MeshType const& ngpMesh,
                                                                 WrappedEntityViewType const& wrappedEntities,
                                                                 int destProc, bool includeGhosts = true)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  auto numWrappedEntities = wrappedEntities.extent(0);
  if (numWrappedEntities == 0) {
    return Kokkos::View<unsigned*>("intersectedIndices", 0);
  }

  Kokkos::View<bool*> entitiesInIntersection("entitiesInIntersection", numWrappedEntities);

  unsigned numEntitiesFound = 0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, numWrappedEntities),
    KOKKOS_LAMBDA(const int i, unsigned& update) {
      auto key = ngpMesh.entity_key(static_cast<Entity>(wrappedEntities(i)));
      auto rank = key.rank();

      bool isRelevantToDestProc = find_entity_key_in_comm_map(ngpMesh, rank, destProc, key,
                                                              includeGhosts).is_local_offset_valid();
      if (isRelevantToDestProc) {
        entitiesInIntersection(i) = true;
        ++update;
      }
    }, Kokkos::Sum<unsigned>(numEntitiesFound)
  );
  Kokkos::fence();

  Kokkos::View<unsigned*> intersectedIndices("intersectedIndices", numEntitiesFound);
  if (numEntitiesFound == 0) { return intersectedIndices; }

  Kokkos::parallel_scan(Kokkos::RangePolicy<ExecSpace>(0, numWrappedEntities),
    KOKKOS_LAMBDA(const size_t i, unsigned& update, const bool isFinal) {
      bool isFoundEntity = entitiesInIntersection(i);
      unsigned pos = update;

      if (isFoundEntity) { ++update; }
      if (isFinal && isFoundEntity) { intersectedIndices(pos) = static_cast<unsigned>(i); }
    }
  );
  Kokkos::fence();

  return intersectedIndices;
}

struct PartChange
{
  PartOrdinal part = InvalidPartOrdinal;
  bool isAdd = false;
};

struct SharedEntityPartChanges
{
  Kokkos::View<EntityKey*, stk::ngp::HostMemSpace> keys;
  Kokkos::View<unsigned*, stk::ngp::HostMemSpace> entityChangeOffsets;
  Kokkos::View<PartChange*, stk::ngp::HostMemSpace> partChanges;
};

struct UnpackedPartChange
{
  EntityKey key;
  PartOrdinal part = InvalidPartOrdinal;
  bool isAdd = false;
  int srcProc = -1;

  bool operator<(UnpackedPartChange const& rhs) const
  {
    if (!(key == rhs.key)) { return key < rhs.key; }
    if (part != rhs.part) { return part < rhs.part; }
    if (isAdd != rhs.isAdd) { return isAdd; }
    return srcProc < rhs.srcProc;
  }
};

template <typename T>
Kokkos::View<T*> copy_vector_to_device_view(std::vector<T> const& vec)
{
  Kokkos::View<T*> deviceView(Kokkos::view_alloc(Kokkos::WithoutInitializing, ""), vec.size());
  Kokkos::View<const T*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> hostView(vec.data(), vec.size());
  Kokkos::deep_copy(deviceView, hostView);
  return deviceView;
}

template <typename OffsetsPerRankType>
KOKKOS_INLINE_FUNCTION
bool is_entity_shared_proc(OffsetsPerRankType const& offsetsPerRank, int numProcs, int proc)
{
  for (int rank = stk::topology::BEGIN_RANK; rank < stk::topology::NUM_RANKS; ++rank) {
    auto const& offsets = offsetsPerRank[rank];
    if (offsets.extent(0) < static_cast<size_t>(numProcs) + 1) { continue; }
    if (offsets(proc+1) > offsets(proc)) { return true; }
  }
  return false;
}

template <typename OffsetsViewType>
Kokkos::View<int*, stk::ngp::HostMemSpace> get_entity_shared_procs(Kokkos::Array<OffsetsViewType, stk::topology::NUM_RANKS> const& offsetsPerRank,
                                                                   int numProcs, int myRank)
{
  using ExecSpace = typename OffsetsViewType::execution_space;

  unsigned numEntitySharedProcs = 0;
  Kokkos::parallel_reduce("count_entity_shared_procs", Kokkos::RangePolicy<ExecSpace>(0, numProcs),
    KOKKOS_LAMBDA(const int proc, unsigned& update) {
      if (proc != myRank && is_entity_shared_proc(offsetsPerRank, numProcs, proc)) { ++update; }
    }, numEntitySharedProcs
  );

  Kokkos::View<int*> entitySharedProcs(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entitySharedProcs"), numEntitySharedProcs);
  if (numEntitySharedProcs > 0) {
    Kokkos::parallel_scan("fill_entity_shared_procs", Kokkos::RangePolicy<ExecSpace>(0, numProcs),
      KOKKOS_LAMBDA(const int proc, unsigned& update, const bool isFinal) {
        if (proc == myRank || !is_entity_shared_proc(offsetsPerRank, numProcs, proc)) { return; }
        if (isFinal) { entitySharedProcs(update) = proc; }
        ++update;
      }
    );
    Kokkos::fence();
  }

  return Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, entitySharedProcs);
}

template <typename PartOrdinalsProxyType, typename PartChangeFunctor>
KOKKOS_INLINE_FUNCTION
void for_each_part_change(Kokkos::pair<const PartOrdinal*, const PartOrdinal*> currentParts, PartOrdinalsProxyType const& newParts,
                          PartChangeFunctor const& partChangeFunctor)
{
  const PartOrdinal* currentPart = currentParts.first;
  const PartOrdinal* newPart = newParts.startPtr;
  const PartOrdinal* const newPartsEnd = newParts.startPtr + newParts.length;
  while (newPart != newPartsEnd || currentPart != currentParts.second) {
    if (currentPart == currentParts.second || (newPart != newPartsEnd && *newPart < *currentPart)) {
      partChangeFunctor(PartChange{*newPart++, true});
    } else if (newPart == newPartsEnd || *currentPart < *newPart) {
      partChangeFunctor(PartChange{*currentPart++, false});
    } else {
      ++newPart;
      ++currentPart;
    }
  }
}

template <typename PartOrdinalsProxyType>
KOKKOS_INLINE_FUNCTION
unsigned get_num_part_changes(Kokkos::pair<const PartOrdinal*, const PartOrdinal*> currentParts, PartOrdinalsProxyType const& newParts)
{
  unsigned numPartChanges = 0;
  for_each_part_change(currentParts, newParts, [&](PartChange const&) { ++numPartChanges; });
  return numPartChanges;
}

template <typename PartOrdinalsProxyType>
KOKKOS_INLINE_FUNCTION
void fill_part_changes(Kokkos::pair<const PartOrdinal*, const PartOrdinal*> currentParts, PartOrdinalsProxyType const& newParts,
                       PartChange* partChanges)
{
  unsigned numPartChanges = 0;
  for_each_part_change(currentParts, newParts, [&](PartChange const& partChange) { partChanges[numPartChanges++] = partChange; });
}

template <typename MeshType, typename NewPartListsType>
SharedEntityPartChanges compute_part_changes_of_shared_entities(MeshType const& ngpMesh,
                                                                NewPartListsType const& newPartLists,
                                                                std::vector<unsigned> const& entityIndices)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  const size_t numSharedEntities = entityIndices.size();

  if (numSharedEntities == 0) { return {}; }

  auto entityIndicesOnDevice = copy_vector_to_device_view(entityIndices);
  auto wrappedEntities = newPartLists.wrappedEntities;
  auto partOrdinalsProxy = newPartLists.partOrdinalsProxy;

  Kokkos::View<EntityKey*> sharedKeys(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sharedKeys"), numSharedEntities);
  Kokkos::View<unsigned*> entityChangeOffsets(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sharedEntityChangeOffsets"), numSharedEntities + 1);

  unsigned numPartChanges = 0;
  Kokkos::parallel_scan("count_part_changes_of_shared_entities", Kokkos::RangePolicy<ExecSpace>(0, numSharedEntities + 1),
    KOKKOS_LAMBDA(const size_t i, unsigned& update, const bool isFinal) {
      if (isFinal) { entityChangeOffsets(i) = update; }
      if (i == numSharedEntities) { return; }

      const unsigned affectedEntityIdx = entityIndicesOnDevice(i);
      const Entity entity = static_cast<Entity>(wrappedEntities(affectedEntityIdx));
      const auto currentParts = ngpMesh.get_bucket(ngpMesh.entity_rank(entity), ngpMesh.device_mesh_index(entity).bucket_id).superset_part_ordinals();
      if (isFinal) { sharedKeys(i) = ngpMesh.entity_key(entity); }
      update += get_num_part_changes(currentParts, partOrdinalsProxy(affectedEntityIdx));
    }, numPartChanges
  );
  Kokkos::fence();

  Kokkos::View<PartChange*> partChanges(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sharedPartChanges"), numPartChanges);
  Kokkos::parallel_for("fill_part_changes_of_shared_entities", Kokkos::RangePolicy<ExecSpace>(0, numSharedEntities),
    KOKKOS_LAMBDA(const size_t i) {
      const unsigned affectedEntityIdx = entityIndicesOnDevice(i);
      const Entity entity = static_cast<Entity>(wrappedEntities(affectedEntityIdx));
      const auto currentParts = ngpMesh.get_bucket(ngpMesh.entity_rank(entity), ngpMesh.device_mesh_index(entity).bucket_id).superset_part_ordinals();
      fill_part_changes(currentParts, partOrdinalsProxy(affectedEntityIdx), partChanges.data() + entityChangeOffsets(i));
    }
  );
  Kokkos::fence();

  SharedEntityPartChanges sharedEntityPartChanges;
  sharedEntityPartChanges.keys = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, sharedKeys);
  sharedEntityPartChanges.entityChangeOffsets = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, entityChangeOffsets);
  sharedEntityPartChanges.partChanges = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, partChanges);
  return sharedEntityPartChanges;
}

// SEND
template <typename MeshType, typename NewPartListsType>
void communicate_entity_part_changes(MeshType& ngpMesh, CommSparse& commSparse,
                                                      Kokkos::View<int*, stk::ngp::HostMemSpace> const& entitySharedProcs,
                                                      NewPartListsType const& localNewPartLists)
{
  Kokkos::Profiling::pushRegion("intersect_wrapped_entities_with_comm_map");
  std::vector<Kokkos::View<unsigned*, stk::ngp::HostMemSpace>> affectedEntityIndicesPerSharedProc(entitySharedProcs.extent(0));
  std::vector<unsigned> affectedEntityIndices;

  for (size_t n = 0; n < entitySharedProcs.extent(0); ++n) {
    auto& indices = affectedEntityIndicesPerSharedProc[n];
    auto entityIndices = intersect_wrapped_entities_with_comm_map(ngpMesh, localNewPartLists.wrappedEntities, entitySharedProcs(n));
    indices = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, entityIndices);
    affectedEntityIndices.insert(affectedEntityIndices.end(), indices.data(), indices.data() + indices.extent(0));
  }

  stk::util::sort_and_unique(affectedEntityIndices);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("compute_part_changes_of_shared_entities");
  SharedEntityPartChanges sharedEntityPartChanges = compute_part_changes_of_shared_entities(ngpMesh, localNewPartLists, affectedEntityIndices);

  std::vector<std::vector<unsigned>> entriesPerEntitySharedProc(entitySharedProcs.extent(0));

  for (size_t n = 0; n < entitySharedProcs.extent(0); ++n) {
    auto& entries = entriesPerEntitySharedProc[n];
    auto& indices = affectedEntityIndicesPerSharedProc[n];
    for (size_t i = 0; i < indices.extent(0); ++i) {
      const unsigned affectedEntityIdx = indices(i);
      const unsigned entry = std::lower_bound(affectedEntityIndices.begin(), affectedEntityIndices.end(), affectedEntityIdx)
                           - affectedEntityIndices.begin();
      if (sharedEntityPartChanges.entityChangeOffsets(entry+1) > sharedEntityPartChanges.entityChangeOffsets(entry)) { entries.push_back(entry); }
    }
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("pack_and_communicate_entity_part_changes");
  stk::pack_and_communicate(commSparse, [&]() {
    for (size_t n = 0; n < entitySharedProcs.extent(0); ++n) {
      auto const& entries = entriesPerEntitySharedProc[n];
      if (entries.empty()) { continue; }

      auto& buf = commSparse.send_buffer(entitySharedProcs(n));
      buf.pack(static_cast<unsigned>(entries.size()));

      for (unsigned entry : entries) {
        const unsigned begin = sharedEntityPartChanges.entityChangeOffsets(entry);
        const unsigned end = sharedEntityPartChanges.entityChangeOffsets(entry+1);

        buf.pack(sharedEntityPartChanges.keys(entry));
        buf.pack(end - begin);
        for (unsigned i = begin; i < end; ++i) {
          buf.pack(sharedEntityPartChanges.partChanges(i).part);
          buf.pack(sharedEntityPartChanges.partChanges(i).isAdd);
        }
      }
    }
  });
  Kokkos::Profiling::popRegion();
}

template <typename MeshType, typename KeyViewType>
Kokkos::View<Entity*> get_entities_from_keys(MeshType const& ngpMesh, KeyViewType const& sortedKeys)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  Kokkos::View<Entity*> localEntities("localEntities", sortedKeys.extent(0));
  if (sortedKeys.extent(0) == 0) { return localEntities; }

  auto& entityKeys = ngpMesh.get_entity_keys();
  const size_t numEntityKeysInMesh = Kokkos::min(entityKeys.extent(0), ngpMesh.get_fast_mesh_indices().extent(0));

  Kokkos::parallel_for("get_entities_from_keys", Kokkos::RangePolicy<ExecSpace>(0, numEntityKeysInMesh),
    KOKKOS_LAMBDA(const size_t localOffset) {
      const EntityKey entityKeyInMesh = entityKeys(localOffset);
      if (!entityKeyInMesh.is_valid()) { return; }

      auto keyIdx = sorted_view_lower_bound(sortedKeys, entityKeyInMesh);
      if (keyIdx < sortedKeys.extent(0) && sortedKeys(keyIdx) == entityKeyInMesh) {
        localEntities(keyIdx) = Entity{static_cast<Entity::entity_value_type>(localOffset)};
      }
    }
  );
  Kokkos::fence();

#ifndef NDEBUG
  unsigned numKeysNotFound = 0;
  Kokkos::parallel_reduce("count_keys_not_found", Kokkos::RangePolicy<ExecSpace>(0, localEntities.extent(0)),
    KOKKOS_LAMBDA(const size_t keyIdx, unsigned& update) {
      if (!localEntities(keyIdx).is_local_offset_valid()) { ++update; }
    }, numKeysNotFound
  );
  STK_ThrowAssertMsg(numKeysNotFound == 0, "Found no local entity for " << numKeysNotFound << " entity keys");
#endif

  return localEntities;
}

inline std::vector<UnpackedPartChange> unpack_received_part_changes(CommSparse& commSparse)
{
  std::vector<UnpackedPartChange> receivedChanges;
  for (int srcProc = 0; srcProc < commSparse.parallel_size(); ++srcProc) {
    if (srcProc == commSparse.parallel_rank()) { continue; }

    stk::CommBuffer& buf = commSparse.recv_buffer(srcProc);
    if (buf.remaining() <= 0) { continue; }

    unsigned numEntities = 0;
    buf.unpack(numEntities);
    for (unsigned e = 0; e < numEntities; ++e) {
      EntityKey key;
      buf.unpack(key);

      unsigned numPartChanges = 0;
      buf.unpack(numPartChanges);
      for (unsigned c = 0; c < numPartChanges; ++c) {
        PartOrdinal part = InvalidPartOrdinal;
        bool isAdd = false;
        buf.unpack(part);
        buf.unpack(isAdd);
        receivedChanges.push_back(UnpackedPartChange{key, part, isAdd, srcProc});
      }
    }
  }
  return receivedChanges;
}

inline bool is_same_entity_part(UnpackedPartChange const& lhs, UnpackedPartChange const& rhs)
{
  return lhs.key == rhs.key && lhs.part == rhs.part;
}

inline void combine_received_part_changes(std::vector<UnpackedPartChange>& receivedChanges)
{
  std::sort(receivedChanges.begin(), receivedChanges.end());

  size_t numKept = 0;
  for (size_t first = 0; first < receivedChanges.size();) {
    size_t last = first + 1;
    while (last < receivedChanges.size() && is_same_entity_part(receivedChanges[last], receivedChanges[first])) { ++last; }

    const size_t keptEnd = receivedChanges[first].isAdd ? first + 1 : last;
    for (size_t c = first; c < keptEnd; ++c) { receivedChanges[numKept++] = receivedChanges[c]; }
    first = last;
  }
  receivedChanges.resize(numKept);
}

// RECV
inline std::vector<UnpackedPartChange> unpack_entity_part_changes(CommSparse& commSparse)
{
  Kokkos::Profiling::pushRegion("unpack_received_part_changes");
  std::vector<UnpackedPartChange> receivedChanges = unpack_received_part_changes(commSparse);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("combine_received_part_changes");
  combine_received_part_changes(receivedChanges);
  Kokkos::Profiling::popRegion();

  return receivedChanges;
}

inline std::vector<EntityKey> get_unique_keys_of_part_changes(std::vector<UnpackedPartChange> const& sortedPartChanges)
{
  std::vector<EntityKey> uniqueKeys;
  for (auto const& partChange : sortedPartChanges) {
    if (uniqueKeys.empty() || !(uniqueKeys.back() == partChange.key)) {
      uniqueKeys.push_back(partChange.key);
    }
  }
  return uniqueKeys;
}

using EntityKeyPartPairs = std::vector<std::pair<EntityKey, PartOrdinal>>;

template <typename MeshType, typename WrappedEntityViewType, typename RmPartsViewType>
EntityKeyPartPairs determine_parts_to_deny_removal(MeshType const& ngpMesh, std::vector<UnpackedPartChange> const& receivedPartChanges,
                                                   Kokkos::View<Entity*> const& receivedEntities,
                                                   WrappedEntityViewType const& wrappedEntities, RmPartsViewType const& rmParts)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using TeamMember = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

  std::vector<unsigned> removalKeyIndices;
  std::vector<PartOrdinal> removalParts;
  std::vector<EntityKey> removalKeys;
  unsigned keyIdx = 0;
  for (size_t i = 0; i < receivedPartChanges.size(); ++i) {
    auto const& partChange = receivedPartChanges[i];
    if (i > 0 && !(receivedPartChanges[i-1].key == partChange.key)) { ++keyIdx; }
    if (partChange.isAdd || (i > 0 && is_same_entity_part(receivedPartChanges[i-1], partChange))) { continue; }

    removalKeyIndices.push_back(keyIdx);
    removalParts.push_back(partChange.part);
    removalKeys.push_back(partChange.key);
  }
  STK_ThrowAssert(receivedPartChanges.empty() || keyIdx + 1 == receivedEntities.extent(0));

  auto numRemovals = removalParts.size();
  if (numRemovals == 0) { return EntityKeyPartPairs{}; }

  auto removalKeyIndicesOnDevice = copy_vector_to_device_view(removalKeyIndices);
  auto removalPartsOnDevice = copy_vector_to_device_view(removalParts);

  Kokkos::View<bool*> shouldDenyPartRemoval("shouldDenyPartRemoval", numRemovals);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for("determine_parts_to_deny_removal", Kokkos::TeamPolicy<ExecSpace>(numRemovals, Kokkos::AUTO),
    KOKKOS_LAMBDA(TeamMember const& team) {
      auto removalIdx = team.league_rank();
      Entity entity = receivedEntities(removalKeyIndicesOnDevice(removalIdx));
      PartOrdinal part = removalPartsOnDevice(removalIdx);
      auto rank = ngpMesh.entity_rank(entity);

      // not ranked part
      if (!deviceBucketRepo.is_inducible_part(part)) { return; }

      auto partRank = deviceBucketRepo.get_part_rank(part);
      // not inducible to current rank
      if (!(partRank > rank)) { return; }

      bool isInRmParts = false;
      for (unsigned rp = 0; rp < rmParts.extent(0); ++rp) {
        if (rmParts(rp) == part) { isInRmParts = true; }
      }

      auto isLosingPart = [&](Entity upperEntity) {
        if (!isInRmParts) { return false; }
        const size_t wrappedEntityIdx = sorted_view_lower_bound(wrappedEntities, upperEntity);
        return wrappedEntityIdx < wrappedEntities.extent(0) && static_cast<Entity>(wrappedEntities(wrappedEntityIdx)) == upperEntity &&
               wrappedEntities(wrappedEntityIdx).isUserInputEntity;
      };

      const auto fastMeshIdx = ngpMesh.device_mesh_index(entity);
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIdx.bucket_id);
      const bool shouldKeepPart = should_keep_induced_part(ngpMesh, team, bucket, fastMeshIdx.bucket_ord,
                                                           part, partRank, isLosingPart);

      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        if (shouldKeepPart) {
          shouldDenyPartRemoval(removalIdx) = true;
        }
      });
    }
  );
  Kokkos::fence();

  auto shouldDenyPartRemovalOnHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, shouldDenyPartRemoval);
  EntityKeyPartPairs deniedPartRemovals;
  for (size_t r = 0; r < numRemovals; ++r) {
    if (shouldDenyPartRemovalOnHost(r)) {
      deniedPartRemovals.emplace_back(removalKeys[r], removalParts[r]);
    }
  }
  return deniedPartRemovals;
}

// Packs and sends each proc its (key, part) pairs
inline void send_entity_key_part_pairs(CommSparse& commSparse, std::vector<EntityKeyPartPairs> const& pairsPerDestProc)
{
  const int numProcs = commSparse.parallel_size();
  const int myRank = commSparse.parallel_rank();

  stk::pack_and_communicate(commSparse, [&]() {
    for (int destProc = 0; destProc < numProcs; ++destProc) {
      auto const& pairs = pairsPerDestProc[destProc];
      if (destProc == myRank || pairs.empty()) { continue; }

      auto& buf = commSparse.send_buffer(destProc);
      buf.pack(static_cast<unsigned>(pairs.size()));
      for (auto const& [key, part] : pairs) {
        buf.pack(key);
        buf.pack(part);
      }
    }
  });
}

// Unpacks the (key, part) pairs received from each proc
inline std::vector<EntityKeyPartPairs> recv_entity_key_part_pairs(CommSparse& commSparse)
{
  const int numProcs = commSparse.parallel_size();
  const int myRank = commSparse.parallel_rank();

  std::vector<EntityKeyPartPairs> pairsPerSrcProc(numProcs);
  for (int srcProc = 0; srcProc < numProcs; ++srcProc) {
    if (srcProc == myRank) { continue; }

    stk::CommBuffer& buf = commSparse.recv_buffer(srcProc);
    if (buf.remaining() <= 0) { continue; }

    unsigned numPairs = 0;
    buf.unpack(numPairs);
    pairsPerSrcProc[srcProc].resize(numPairs);
    for (auto& [key, part] : pairsPerSrcProc[srcProc]) {
      buf.unpack(key);
      buf.unpack(part);
    }
  }
  return pairsPerSrcProc;
}


inline std::vector<EntityKeyPartPairs> collect_denied_part_removals_per_src_proc(std::vector<UnpackedPartChange> const& receivedPartChanges,
                                                                                EntityKeyPartPairs const& deniedPartRemovals, int numProcs)
{
  std::vector<EntityKeyPartPairs> deniedPartRemovalsPerSrcProc(numProcs);
  if (deniedPartRemovals.empty()) { return deniedPartRemovalsPerSrcProc; }

  for (auto const& partChange : receivedPartChanges) {
    if (!partChange.isAdd &&
        std::binary_search(deniedPartRemovals.begin(), deniedPartRemovals.end(), std::make_pair(partChange.key, partChange.part))) {
      deniedPartRemovalsPerSrcProc[partChange.srcProc].emplace_back(partChange.key, partChange.part);
    }
  }
  return deniedPartRemovalsPerSrcProc;
}

inline void drop_denied_part_removals(std::vector<UnpackedPartChange>& receivedPartChanges, EntityKeyPartPairs const& deniedPartRemovals)
{
  if (deniedPartRemovals.empty()) { return; }

  std::erase_if(receivedPartChanges, [&](UnpackedPartChange const& partChange) {
    return !partChange.isAdd &&
           std::binary_search(deniedPartRemovals.begin(), deniedPartRemovals.end(), std::make_pair(partChange.key, partChange.part));
  });
}

template <typename MeshType, typename WrappedEntityViewType, typename RmPartsViewType>
void send_denied_part_removals(CommSparse& commSparse, MeshType const& ngpMesh, std::vector<UnpackedPartChange>& receivedPartChanges,
                               Kokkos::View<Entity*> const& receivedEntities,
                               WrappedEntityViewType const& wrappedEntities, RmPartsViewType const& rmParts)
{
  Kokkos::Profiling::pushRegion("determine_parts_to_deny_removal");
  auto localDeniedPartRemovals = determine_parts_to_deny_removal(ngpMesh, receivedPartChanges, receivedEntities, wrappedEntities, rmParts);
  auto deniedPartRemovalsPerSrcProc = collect_denied_part_removals_per_src_proc(receivedPartChanges, localDeniedPartRemovals,
                                                                               commSparse.parallel_size());
  drop_denied_part_removals(receivedPartChanges, localDeniedPartRemovals);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("pack_and_communicate_denied_part_removals");
  send_entity_key_part_pairs(commSparse, deniedPartRemovalsPerSrcProc);
  Kokkos::Profiling::popRegion();
}

inline EntityKeyPartPairs recv_denied_part_removals(CommSparse& commSparse)
{
  EntityKeyPartPairs deniedKeyPartPairs;
  for (auto const& received : recv_entity_key_part_pairs(commSparse)) {
    deniedKeyPartPairs.insert(deniedKeyPartPairs.end(), received.begin(), received.end());
  }
  stk::util::sort_and_unique(deniedKeyPartPairs);
  return deniedKeyPartPairs;
}

// Puts the parts that were denied from removing by at least one other proc back into local to-be-removed part list
template <typename MeshType, typename WrappedEntityViewType, typename RmPartsViewType,
          typename PartOrdinalsProxyViewType, typename DeniedPartRemovalViewType>
void restore_denied_part_removals([[maybe_unused]] MeshType const& ngpMesh,
                                  WrappedEntityViewType const& wrappedEntities,
                                  [[maybe_unused]] RmPartsViewType const& rmParts,
                                  [[maybe_unused]] unsigned maxPartsPerEntity,
                                  PartOrdinalsProxyViewType const& partOrdinalsProxyView,
                                  DeniedPartRemovalViewType const& deniedPartRemovalView)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  const size_t numDenied = deniedPartRemovalView.extent(0);
  if (numDenied == 0) { return; }

#ifndef NDEBUG
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  Kokkos::parallel_for("check_denied_part_removals", Kokkos::RangePolicy<ExecSpace>(0, numDenied),
    KOKKOS_LAMBDA(const size_t deniedIdx) {
      const Entity entity = deniedPartRemovalView(deniedIdx).entity;
      const Ordinal part = deniedPartRemovalView(deniedIdx).part;

      const size_t entityIdx = sorted_view_lower_bound(wrappedEntities, entity);
      STK_NGP_ThrowAssert(entityIdx < wrappedEntities.extent(0) && static_cast<Entity>(wrappedEntities(entityIdx)) == entity);

      const auto rank = ngpMesh.entity_rank(entity);
      const auto currentParts = ngpMesh.get_bucket(rank, ngpMesh.device_mesh_index(entity).bucket_id).superset_part_ordinals();
      bool isInRmParts = false;
      for (unsigned rp = 0; rp < rmParts.extent(0); ++rp) {
        if (rmParts(rp) == part) { isInRmParts = true; }
      }
      bool isCurrentPart = false;
      for (const Ordinal* currentPart = currentParts.first; currentPart != currentParts.second; ++currentPart) {
        if (*currentPart == part) { isCurrentPart = true; }
      }
      STK_NGP_ThrowAssert(deviceBucketRepo.is_inducible_part(part));
      STK_NGP_ThrowAssert(deviceBucketRepo.get_part_rank(part) > rank);
      STK_NGP_ThrowAssert(isInRmParts);
      STK_NGP_ThrowAssert(isCurrentPart);
    }
  );
  Kokkos::fence();
#endif

  Kokkos::parallel_for("restore_denied_part_removals", Kokkos::RangePolicy<ExecSpace>(0, numDenied),
    KOKKOS_LAMBDA(const size_t firstDeniedIdx) {
      const Entity entity = deniedPartRemovalView(firstDeniedIdx).entity;
      if (firstDeniedIdx > 0 && deniedPartRemovalView(firstDeniedIdx-1).entity == entity) { return; }

      const size_t entityIdx = sorted_view_lower_bound(wrappedEntities, entity);
      Ordinal* newParts = partOrdinalsProxyView(entityIdx).startPtr;
      unsigned numNewParts = partOrdinalsProxyView(entityIdx).length;

      for (size_t deniedIdx = firstDeniedIdx; deniedIdx < numDenied && deniedPartRemovalView(deniedIdx).entity == entity; ++deniedIdx) {
        const Ordinal part = deniedPartRemovalView(deniedIdx).part;

        unsigned insertIdx = 0;
        while (insertIdx < numNewParts && newParts[insertIdx] < part) { ++insertIdx; }
        if (insertIdx < numNewParts && newParts[insertIdx] == part) { continue; }

        STK_NGP_ThrowAssert(numNewParts < maxPartsPerEntity);
        for (unsigned i = numNewParts; i > insertIdx; --i) {
          newParts[i] = newParts[i-1];
        }
        newParts[insertIdx] = part;
        ++numNewParts;
      }

      partOrdinalsProxyView(entityIdx).length = numNewParts;
    }
  );
  Kokkos::fence();
}

// Forwards the denials of this proc's removals also to the other procs holding those entities
template <typename MeshType>
void forward_denied_part_removals(MeshType& ngpMesh, std::vector<UnpackedPartChange>& receivedChanges,
                                  Kokkos::View<int*, stk::ngp::HostMemSpace> const& entitySharedProcs,
                                  EntityKeyPartPairs const& deniedKeyPartPairs, std::vector<unsigned> const& deniedKeyOffsets,
                                  Kokkos::View<Entity*> const& deniedEntities)
{
  const ParallelMachine comm = ngpMesh.get_bulk_on_host().parallel();
  const int numProcs = stk::parallel_machine_size(comm);

  std::vector<EntityKeyPartPairs> deniedPartRemovalsPerProc(numProcs);

  if (deniedEntities.extent(0) > 0) {
    for (size_t n = 0; n < entitySharedProcs.extent(0); ++n) {
      int proc = entitySharedProcs(n);
      auto deniedIndicesOfProc = intersect_wrapped_entities_with_comm_map(ngpMesh, deniedEntities, proc);
      auto deniedIndices = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, deniedIndicesOfProc);
      for (size_t i = 0; i < deniedIndices.extent(0); ++i) {
        const unsigned deniedIdx = deniedIndices(i);
        deniedPartRemovalsPerProc[proc].insert(deniedPartRemovalsPerProc[proc].end(),
                                               deniedKeyPartPairs.begin() + deniedKeyOffsets[deniedIdx],
                                               deniedKeyPartPairs.begin() + deniedKeyOffsets[deniedIdx + 1]);
      }
    }
  }

  stk::CommSparse commSparse(comm);
  send_entity_key_part_pairs(commSparse, deniedPartRemovalsPerProc);
  drop_denied_part_removals(receivedChanges, recv_denied_part_removals(commSparse));
}

// Start of each key's run of pairs in sorted (key, part) pairs, plus the end
inline std::vector<unsigned> get_key_offsets(EntityKeyPartPairs const& sortedKeyPartPairs)
{
  std::vector<unsigned> keyOffsets;
  for (size_t i = 0; i < sortedKeyPartPairs.size(); ++i) {
    if (i == 0 || !(sortedKeyPartPairs[i].first == sortedKeyPartPairs[i-1].first)) { keyOffsets.push_back(i); }
  }
  keyOffsets.push_back(sortedKeyPartPairs.size());
  return keyOffsets;
}

inline Kokkos::View<EntityKey*, stk::ngp::HostMemSpace> get_unique_keys(EntityKeyPartPairs const& sortedKeyPartPairs,
                                                                       std::vector<unsigned> const& keyOffsets)
{
  Kokkos::View<EntityKey*, stk::ngp::HostMemSpace> uniqueKeys(Kokkos::view_alloc(Kokkos::WithoutInitializing, "uniqueKeys"),
                                                              keyOffsets.size() - 1);
  for (size_t keyIdx = 0; keyIdx < uniqueKeys.extent(0); ++keyIdx) {
    uniqueKeys(keyIdx) = sortedKeyPartPairs[keyOffsets[keyIdx]].first;
  }
  return uniqueKeys;
}

// Device view of the denied (entity, part) removals
template <typename DeniedViewType>
DeniedViewType make_denied_part_removal_view(EntityKeyPartPairs const& deniedKeyPartPairs, std::vector<unsigned> const& deniedKeyOffsets,
                                             Kokkos::View<Entity*> const& deniedEntities)
{
  auto deniedEntitiesOnHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, deniedEntities);

  Kokkos::View<DeniedPartRemoval*, stk::ngp::HostMemSpace> deniedPartRemovalsOnHost(Kokkos::view_alloc(Kokkos::WithoutInitializing, "deniedPartRemovalsOnHost"), deniedKeyPartPairs.size());

  for (size_t deniedIdx = 0; deniedIdx < deniedEntitiesOnHost.extent(0); ++deniedIdx) {
    for (unsigned pairIdx = deniedKeyOffsets[deniedIdx]; pairIdx < deniedKeyOffsets[deniedIdx + 1]; ++pairIdx) {
      deniedPartRemovalsOnHost(pairIdx) = DeniedPartRemoval{deniedEntitiesOnHost(deniedIdx), deniedKeyPartPairs[pairIdx].second};
    }
  }
  std::sort(deniedPartRemovalsOnHost.data(), deniedPartRemovalsOnHost.data() + deniedPartRemovalsOnHost.extent(0));

  return Kokkos::create_mirror_view_and_copy(typename DeniedViewType::memory_space{}, deniedPartRemovalsOnHost);
}

// Deny Part Removals
template <typename MeshType, typename WrappedEntityViewType, typename RmPartsViewType>
Kokkos::View<DeniedPartRemoval*, typename MeshType::ngp_mem_space>
communicate_denied_part_removals(MeshType& ngpMesh, std::vector<UnpackedPartChange>& receivedPartChanges,
                                 Kokkos::View<Entity*> const& receivedEntities,
                                 Kokkos::View<int*, stk::ngp::HostMemSpace> const& entitySharedProcs,
                                 WrappedEntityViewType const& wrappedEntities, RmPartsViewType const& rmParts)
{
  using DeniedViewType = Kokkos::View<DeniedPartRemoval*, typename MeshType::ngp_mem_space>;

  stk::CommSparse commSparse(ngpMesh.get_bulk_on_host().parallel());
  send_denied_part_removals(commSparse, ngpMesh, receivedPartChanges, receivedEntities, wrappedEntities, rmParts);

  Kokkos::Profiling::pushRegion("unpack_denied_part_removals");
  auto deniedKeyPartPairs = recv_denied_part_removals(commSparse);
  auto deniedKeyOffsets = get_key_offsets(deniedKeyPartPairs);
  auto deniedKeys = get_unique_keys(deniedKeyPartPairs, deniedKeyOffsets);
  auto deniedEntities = get_entities_from_keys(ngpMesh, Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, deniedKeys));
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("forward_denied_part_removals");
  forward_denied_part_removals(ngpMesh, receivedPartChanges, entitySharedProcs, deniedKeyPartPairs, deniedKeyOffsets, deniedEntities);
  Kokkos::Profiling::popRegion();

  return make_denied_part_removal_view<DeniedViewType>(deniedKeyPartPairs, deniedKeyOffsets, deniedEntities);
}

// New part list of each received entity. own new list or current parts + received changes
template <typename MeshType, typename WrappedEntityViewType, typename OffsetViewType, typename ChangeViewType,
          typename LocalIndexViewType, typename LocalPartOrdinalsProxyViewType,
          typename NewPartsViewType, typename PartOrdinalsProxyViewType>
void set_new_part_list_per_entity_from_part_changes(MeshType const& ngpMesh,
                                                    WrappedEntityViewType const& entities,
                                                    OffsetViewType const& entityChangeOffsets,
                                                    ChangeViewType const& partChanges,
                                                    LocalIndexViewType const& localIndexOfReceivedEntity,
                                                    LocalPartOrdinalsProxyViewType const& localPartOrdinalsProxy,
                                                    unsigned maxPartsPerEntity,
                                                    NewPartsViewType const& newPartsPerEntity,
                                                    PartOrdinalsProxyViewType const& partOrdinalsProxyView)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  Kokkos::parallel_for("set_new_part_lists_from_part_changes", Kokkos::RangePolicy<ExecSpace>(0, entities.extent(0)),
    KOKKOS_LAMBDA(const size_t idx) {
      const Entity entity = static_cast<Entity>(entities(idx));
      const auto rank = ngpMesh.entity_rank(entity);

      const Ordinal* currentPart = nullptr;
      const Ordinal* currentPartsEnd = nullptr;
      const unsigned localIdx = localIndexOfReceivedEntity(idx);
      if (localIdx != INVALID_INDEX) {
        currentPart = localPartOrdinalsProxy(localIdx).startPtr;
        currentPartsEnd = currentPart + localPartOrdinalsProxy(localIdx).length;
      } else {
        const auto currentParts = ngpMesh.get_bucket(rank, ngpMesh.device_mesh_index(entity).bucket_id).superset_part_ordinals();
        currentPart = currentParts.first;
        currentPartsEnd = currentParts.second;
      }

      Ordinal* dest = &newPartsPerEntity(idx * maxPartsPerEntity);
      unsigned numNewParts = 0;

      unsigned changeIdx = entityChangeOffsets(idx);
      const unsigned changesEnd = entityChangeOffsets(idx + 1);
      while (currentPart != currentPartsEnd || changeIdx < changesEnd) {
        if (changeIdx == changesEnd || (currentPart != currentPartsEnd && *currentPart < partChanges(changeIdx).part)) {
          dest[numNewParts++] = *currentPart++;
        } else if (currentPart == currentPartsEnd || partChanges(changeIdx).part < *currentPart) {
          if (partChanges(changeIdx).isAdd) { dest[numNewParts++] = partChanges(changeIdx).part; }
          ++changeIdx;
        } else {
          if (partChanges(changeIdx).isAdd) { dest[numNewParts++] = *currentPart; }
          ++currentPart;
          ++changeIdx;
        }
      }
      STK_NGP_ThrowAssert(numNewParts <= maxPartsPerEntity);

      partOrdinalsProxyView(idx) = typename PartOrdinalsProxyViewType::value_type(rank, dest, numNewParts);
    }
  );
  Kokkos::fence();
}

// Groups the sorted received part changes by entity (index into allReceivedKeys, CSR offsets, changes); removes duplicates
inline unsigned group_received_part_changes_by_entity(std::vector<UnpackedPartChange> const& receivedPartChanges,
                                                      std::vector<EntityKey> const& allReceivedKeys,
                                                      std::vector<unsigned>& receivedKeyIndices,
                                                      std::vector<unsigned>& entityChangeOffsets,
                                                      std::vector<PartChange>& partChanges)
{
  unsigned maxNumAddedParts = 0;
  unsigned numAddedParts = 0;
  size_t keyIdx = 0;
  entityChangeOffsets.assign(1, 0);
  for (auto const& partChange : receivedPartChanges) {
    const bool isNewKey = receivedKeyIndices.empty() || !(allReceivedKeys[receivedKeyIndices.back()] == partChange.key);
    if (isNewKey) {
      if (!receivedKeyIndices.empty()) { entityChangeOffsets.push_back(partChanges.size()); }
      while (!(allReceivedKeys[keyIdx] == partChange.key)) { ++keyIdx; }
      receivedKeyIndices.push_back(keyIdx);
      numAddedParts = 0;
    } else if (partChanges.back().part == partChange.part) {
      continue;
    }

    partChanges.push_back(PartChange{partChange.part, partChange.isAdd});
    if (partChange.isAdd) { maxNumAddedParts = std::max(maxNumAddedParts, ++numAddedParts); }
  }
  entityChangeOffsets.push_back(partChanges.size());
  return maxNumAddedParts;
}

template <typename ExecSpace>
Kokkos::View<Entity*> select_entities_at_indices(Kokkos::View<Entity*> const& entities, std::vector<unsigned> const& indicesOnHost)
{
  auto indices = copy_vector_to_device_view(indicesOnHost);
  Kokkos::View<Entity*> selectedEntities(Kokkos::view_alloc(Kokkos::WithoutInitializing, "selectedEntities"), indices.extent(0));
  Kokkos::parallel_for("select_entities_at_indices", Kokkos::RangePolicy<ExecSpace>(0, indices.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      selectedEntities(i) = entities(indices(i));
    }
  );
  Kokkos::fence();
  return selectedEntities;
}

// Index of each received entity among the local entities, INVALID_INDEX if it has no local entry
template <typename ExecSpace, typename WrappedEntityViewType>
void find_local_indices_of_received_entities(WrappedEntityViewType const& localEntities, Kokkos::View<Entity*> const& receivedEntities,
                                             Kokkos::View<unsigned*> const& localIndexOfReceivedEntity)
{
  Kokkos::parallel_for("find_local_indices_of_received_entities", Kokkos::RangePolicy<ExecSpace>(0, receivedEntities.extent(0)),
    KOKKOS_LAMBDA(const size_t receivedIdx) {
      const Entity entity = receivedEntities(receivedIdx);
      const size_t localIdx = sorted_view_lower_bound(localEntities, entity);
      const bool isLocal = localIdx < localEntities.extent(0) && static_cast<Entity>(localEntities(localIdx)) == entity;
      localIndexOfReceivedEntity(receivedIdx) = isLocal ? static_cast<unsigned>(localIdx) : INVALID_INDEX;
    }
  );
  Kokkos::fence();
}

template <typename ExecSpace, typename PartOrdinalsProxyViewType>
void replace_local_part_lists_with_merged_part_lists(PartOrdinalsProxyViewType const& localPartOrdinalsProxy,
                                                     Kokkos::View<unsigned*> const& localIndexOfReceivedEntity,
                                                     PartOrdinalsProxyViewType const& receivedPartOrdinalsProxy)
{
  Kokkos::parallel_for("replace_local_part_lists_with_merged_part_lists", Kokkos::RangePolicy<ExecSpace>(0, localIndexOfReceivedEntity.extent(0)),
    KOKKOS_LAMBDA(const size_t receivedIdx) {
      const unsigned localIdx = localIndexOfReceivedEntity(receivedIdx);
      if (localIdx != INVALID_INDEX) { localPartOrdinalsProxy(localIdx) = receivedPartOrdinalsProxy(receivedIdx); }
    }
  );
  Kokkos::fence();
}

// Appends the received entities to the currently running adjusted entity/new-parts list
template <typename ExecSpace, typename NewPartListsType, typename PartOrdinalsProxyViewType>
void append_received_only_entities(NewPartListsType& localNewPartLists, Kokkos::View<unsigned*> const& localIndexOfReceivedEntity,
                                   Kokkos::View<Entity*> const& receivedEntities, PartOrdinalsProxyViewType const& receivedPartOrdinalsProxy)
{
  const auto localIndexOfReceivedEntityOnHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, localIndexOfReceivedEntity);
  std::vector<unsigned> receivedOnlyIndicesOnHost;
  for (size_t receivedIdx = 0; receivedIdx < localIndexOfReceivedEntityOnHost.extent(0); ++receivedIdx) {
    if (localIndexOfReceivedEntityOnHost(receivedIdx) == INVALID_INDEX) { receivedOnlyIndicesOnHost.push_back(receivedIdx); }
  }
  if (receivedOnlyIndicesOnHost.empty()) { return; }

  const size_t numLocalEntities = localNewPartLists.wrappedEntities.extent(0);
  const size_t numReceivedOnly = receivedOnlyIndicesOnHost.size();
  Kokkos::resize(localNewPartLists.wrappedEntities, numLocalEntities + numReceivedOnly);
  Kokkos::resize(localNewPartLists.partOrdinalsProxy, numLocalEntities + numReceivedOnly);

  auto receivedOnlyIndices = copy_vector_to_device_view(receivedOnlyIndicesOnHost);
  auto wrappedEntities = localNewPartLists.wrappedEntities;
  auto partOrdinalsProxy = localNewPartLists.partOrdinalsProxy;
  Kokkos::parallel_for("append_received_only_entities", Kokkos::RangePolicy<ExecSpace>(0, numReceivedOnly),
    KOKKOS_LAMBDA(const size_t i) {
      const unsigned receivedIdx = receivedOnlyIndices(i);
      wrappedEntities(numLocalEntities + i) = EntityWrapper<>(receivedEntities(receivedIdx));
      partOrdinalsProxy(numLocalEntities + i) = receivedPartOrdinalsProxy(receivedIdx);
    }
  );
  Kokkos::fence();
}

template <typename MeshType, typename NewPartListsType>
Kokkos::View<Ordinal*> merge_received_part_changes_into_part_lists(MeshType const& ngpMesh, NewPartListsType& localNewPartLists,
                                                                    std::vector<UnpackedPartChange> const& receivedPartChanges,
                                                                    std::vector<EntityKey> const& allReceivedKeys,
                                                                    Kokkos::View<Entity*> const& allReceivedEntities)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using PartOrdinalsProxyViewType = std::remove_cvref_t<decltype(localNewPartLists.partOrdinalsProxy)>;

  std::vector<unsigned> receivedKeyIndicesOnHost;
  std::vector<unsigned> entityChangeOffsetsOnHost;
  std::vector<PartChange> partChangesOnHost;
  const unsigned maxNumAddedParts = group_received_part_changes_by_entity(receivedPartChanges, allReceivedKeys, receivedKeyIndicesOnHost,
                                                                          entityChangeOffsetsOnHost, partChangesOnHost);

  const size_t numReceived = receivedKeyIndicesOnHost.size();
  auto receivedEntities = select_entities_at_indices<ExecSpace>(allReceivedEntities, receivedKeyIndicesOnHost);

  Kokkos::Profiling::pushRegion("find_local_indices_of_received_entities");
  Kokkos::View<unsigned*> localIndexOfReceivedEntity(Kokkos::view_alloc(Kokkos::WithoutInitializing, "localIndexOfReceivedEntity"), numReceived);
  find_local_indices_of_received_entities<ExecSpace>(localNewPartLists.wrappedEntities, receivedEntities, localIndexOfReceivedEntity);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("set_new_part_lists_from_part_changes");
  const unsigned maxNewNumPartsPerEntity = std::max(get_max_num_parts_per_entity(ngpMesh, receivedEntities), localNewPartLists.maxNewNumPartsPerEntity) + maxNumAddedParts;

  Kokkos::View<Ordinal*> receivedNewPartOrdinals("receivedNewPartOrdinals", numReceived * maxNewNumPartsPerEntity);
  PartOrdinalsProxyViewType receivedPartOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "receivedPartOrdinalsProxy"), numReceived);

  set_new_part_list_per_entity_from_part_changes(ngpMesh, receivedEntities, copy_vector_to_device_view(entityChangeOffsetsOnHost),
                                                 copy_vector_to_device_view(partChangesOnHost), localIndexOfReceivedEntity,
                                                 localNewPartLists.partOrdinalsProxy, maxNewNumPartsPerEntity,
                                                 receivedNewPartOrdinals, receivedPartOrdinalsProxy);
  replace_local_part_lists_with_merged_part_lists<ExecSpace>(localNewPartLists.partOrdinalsProxy, localIndexOfReceivedEntity, receivedPartOrdinalsProxy);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("append_received_only_entities");
  append_received_only_entities<ExecSpace>(localNewPartLists, localIndexOfReceivedEntity, receivedEntities, receivedPartOrdinalsProxy);
  Kokkos::Profiling::popRegion();

  return receivedNewPartOrdinals;
}

template <typename MeshType, typename NewPartListsType, typename... RemovePartParams>
Kokkos::View<Ordinal*> communicate_and_merge_part_changes(MeshType& ngpMesh, NewPartListsType& localNewPartLists,
                                                          const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
{
  Kokkos::View<Ordinal*> receivedNewPartOrdinals;
  auto comm = ngpMesh.get_bulk_on_host().parallel();
  stk::CommSparse commSparse(comm);

  Kokkos::Profiling::pushRegion("get_entity_shared_procs");
  const auto entitySharedProcs = get_entity_shared_procs(ngpMesh.volatile_fast_shared_comm_map_offsets(), commSparse.parallel_size(),
                                                         commSparse.parallel_rank());
  Kokkos::Profiling::popRegion();

  // Send part changes of owned + sharing entities to shared procs
  Kokkos::Profiling::pushRegion("communicate_entity_part_changes");
  communicate_entity_part_changes(ngpMesh, commSparse, entitySharedProcs, localNewPartLists);
  Kokkos::Profiling::popRegion();

  // Receive and combine other procs' part changes of entities this proc holds
  Kokkos::Profiling::pushRegion("unpack_entity_part_changes");
  auto receivedPartChanges = unpack_entity_part_changes(commSparse);
  const auto receivedKeys = get_unique_keys_of_part_changes(receivedPartChanges);
  const auto receivedEntities = get_entities_from_keys(ngpMesh, copy_vector_to_device_view(receivedKeys));
  Kokkos::Profiling::popRegion();

  // Determine parts to deny removal (from sharing procs)
  Kokkos::Profiling::pushRegion("all_reduce_has_ranked_remove_part");
  bool hasRankedPartToRemove = has_ranked_remove_part(ngpMesh.get_device_bucket_repository(), removePartOrdinals);
  bool anyProcHasRankedPartToRemove = stk::is_true_on_any_proc(comm, hasRankedPartToRemove);
  Kokkos::Profiling::popRegion();

  if (anyProcHasRankedPartToRemove) {
    Kokkos::Profiling::pushRegion("communicate_denied_part_removals");
    auto deniedPartRemovalView = communicate_denied_part_removals(ngpMesh, receivedPartChanges, receivedEntities, entitySharedProcs,
                                                                        localNewPartLists.wrappedEntities, removePartOrdinals);
    Kokkos::Profiling::popRegion();

    if (deniedPartRemovalView.extent(0) > 0) {
      // Put denied parts back into the local new part lists
      Kokkos::Profiling::pushRegion("restore_denied_part_removals");
      restore_denied_part_removals(ngpMesh, localNewPartLists.wrappedEntities, removePartOrdinals,
                                         localNewPartLists.maxNewNumPartsPerEntity, localNewPartLists.partOrdinalsProxy,
                                         deniedPartRemovalView);
      Kokkos::Profiling::popRegion();
    }
  }

  if (!receivedPartChanges.empty()) {
    // Merge received parts to change and local parts to change
    Kokkos::Profiling::pushRegion("merge_received_part_changes");
    receivedNewPartOrdinals = merge_received_part_changes_into_part_lists(ngpMesh, localNewPartLists, receivedPartChanges,
                                                                             receivedKeys, receivedEntities);
    Kokkos::Profiling::popRegion();
  }

  return receivedNewPartOrdinals;
}

} // namespace stk::mesh::impl

#endif
