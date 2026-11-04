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

#ifndef STK_MESH_NGPMESHIMPL_HPP
#define STK_MESH_NGPMESHIMPL_HPP

#include "stk_util/stk_config.h"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/baseImpl/ViewVector.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include "Kokkos_Sort.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Macros.hpp"
#include <iostream>

namespace stk {
namespace mesh {
namespace impl {

constexpr unsigned INVALID_INDEX = std::numeric_limits<unsigned>::max();

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

struct DevicePartOrdinalLess
{
  KOKKOS_DEFAULTED_FUNCTION
  DevicePartOrdinalLess() = default;

  template <typename PartOrdinal>
  KOKKOS_INLINE_FUNCTION
  bool operator()(PartOrdinal const& lhs, PartOrdinal const& rhs) {
    if (lhs.extent(0) != rhs.extent(0)) {
      return lhs.extent(0) < rhs.extent(0);
    }

    for (unsigned i = 0; i < lhs.extent(0); ++i) {
      if (lhs(i) != rhs(i)) {
        return lhs(i) < rhs(i);
      }
    }
    return false;
  }

  bool operator()(const std::vector<stk::mesh::PartOrdinal>& lhs, const std::vector<stk::mesh::PartOrdinal>& rhs)
  {
    using VectorWrapper = Kokkos::View<const stk::mesh::PartOrdinal*, stk::ngp::HostSpace::mem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    VectorWrapper lhsWrapper(lhs.data(), lhs.size());
    VectorWrapper rhsWrapper(rhs.data(), rhs.size());

    return (*this)(lhsWrapper, rhsWrapper);
  }
};

template <typename T = Entity>
struct EntityWrapper
{
  KOKKOS_FUNCTION
  EntityWrapper()
    : value(T{}),
      isUserInputEntity(false),
      isForPartInduction(false)
  {}

  KOKKOS_FUNCTION
  EntityWrapper(T e, bool u = false, bool p = false)
    : value(e),
      isUserInputEntity(u),
      isForPartInduction(p)
  {}

  KOKKOS_INLINE_FUNCTION
  operator T() const {
    return value;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(EntityWrapper other) const {
    return (value == other.value);
  }

  KOKKOS_INLINE_FUNCTION
  bool operator<(EntityWrapper other) const {
    if (value == other.value) {
      return isForPartInduction == false;
    } else {
      return value < other.value;
    }
  }

  T value;
  bool isUserInputEntity;
  bool isForPartInduction;
};

struct DeniedPartRemoval
{
  Entity entity;
  PartOrdinal part;

  KOKKOS_INLINE_FUNCTION
  bool operator<(DeniedPartRemoval const& rhs) const {
    if (entity.local_offset() != rhs.entity.local_offset()) {
      return entity.local_offset() < rhs.entity.local_offset();
    }
    return part < rhs.part;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(DeniedPartRemoval const& rhs) const {
    return (entity == rhs.entity) && (part == rhs.part);
  }
};

struct NumNewBucketsToAddPerPartition {
  EntityRank rank = topology::INVALID_RANK;
  unsigned partitionId = INVALID_INDEX;
  unsigned numBucketsToAdd = 0;

  KOKKOS_INLINE_FUNCTION
  bool operator==(NumNewBucketsToAddPerPartition const& rhs) const
  {
    return (rank == rhs.rank) && (partitionId == rhs.partitionId) && (numBucketsToAdd == rhs.numBucketsToAdd);
  }
};

struct GrowLastBucketInPartition {
  EntityRank rank = topology::INVALID_RANK;
  unsigned partitionId = INVALID_INDEX;
  bool growLastBucket = false;

  KOKKOS_INLINE_FUNCTION
  bool operator==(GrowLastBucketInPartition const& rhs) const
  {
    return (rank == rhs.rank) && (partitionId == rhs.partitionId) && (growLastBucket == rhs.growLastBucket);
  }
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const EntityWrapper<T>& entityWrapper)
{
  os << static_cast<T>(entityWrapper);
  return os;
}

template<typename MESH_TYPE, typename EntityViewType>
unsigned get_max_num_parts_per_entity(const MESH_TYPE& ngpMesh, const EntityViewType& entities)
{
  using BucketType = typename MESH_TYPE::BucketType;

  unsigned max = 0;

  Kokkos::parallel_reduce( "MaxReduce", entities.size(),
    KOKKOS_LAMBDA (const int& i, unsigned& lmax) {
      const stk::mesh::EntityRank rank = ngpMesh.entity_rank(entities(i));
      const stk::mesh::FastMeshIndex entityIdx = ngpMesh.device_mesh_index(entities(i));
      const BucketType& bucket = ngpMesh.get_bucket(rank, entityIdx.bucket_id);
      auto partOrdinalsPair = bucket.superset_part_ordinals();
      auto numParts = partOrdinalsPair.second - partOrdinalsPair.first;
      lmax = (numParts > lmax) ? numParts : lmax;
    }, Kokkos::Max<unsigned>(max)
  );
  Kokkos::fence();

  return max;
}

template <typename EntityViewType>
Kokkos::View<EntityWrapper<>*> wrap_entities(EntityViewType entities, bool isUserInputEntity = false, bool isForPartInduction = false)
{
  Kokkos::View<EntityWrapper<>*> v("", entities.extent(0));
  
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int i) {
      auto entity = entities(i);
      auto wrappedEntity = EntityWrapper<>{entity, isUserInputEntity, isForPartInduction};
      v(i) = wrappedEntity;
    }
  );
  return v;
}

template <typename EntityView, typename EntityKeyWrapperView>
Kokkos::View<EntityWrapper<>*> convert_to_wrapped_entities(EntityView entities, EntityKeyWrapperView wrappedEntityKeys)
{
  Kokkos::View<EntityWrapper<>*> v("", entities.extent(0));

  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int i) {
      auto entity = entities(i);
      auto wrappedEntityKey = wrappedEntityKeys(i);
      auto isUserInputEntity = wrappedEntityKey.isUserInputEntity;
      auto isForPartInduction = wrappedEntityKey.isForPartInduction;

      v(i) = EntityWrapper<>(entity, isUserInputEntity, isForPartInduction);
    }
  );
  return v;
}

struct LocallyAvailEntityKeys
{
  Kokkos::View<EntityWrapper<EntityKey>*> keys;
  Kokkos::View<Entity*> entities;
};

template <typename MeshType, typename EntityKeyWrapperView>
LocallyAvailEntityKeys get_locally_avail_entity_keys(MeshType const& ngpMesh, EntityKeyWrapperView wrappedEntityKeys)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using TeamMember = typename TeamPolicy::member_type;

  auto numEntityKeys = wrappedEntityKeys.extent(0);
  Kokkos::View<Entity*> tempEntities("tempEntities", numEntityKeys);

  auto& entityKeys = ngpMesh.get_entity_keys();
  auto numEntityKeysInMesh = Kokkos::min(entityKeys.extent(0), ngpMesh.get_fast_mesh_indices().extent(0));

  unsigned numValidEntityKeys = 0;
  Kokkos::parallel_reduce(TeamPolicy(numEntityKeysInMesh, Kokkos::AUTO),
    KOKKOS_LAMBDA(TeamMember const& team, unsigned& update) {
      auto i = team.league_rank();
      auto entityKeyInMesh = entityKeys(i);
      if (!entityKeyInMesh.is_valid()) { return; }

      bool found = false;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, numEntityKeys),
        [&](const int j) {
          if (found) { return; }
          if (wrappedEntityKeys(j).value == entityKeyInMesh) {
            tempEntities(j) = Entity{static_cast<unsigned>(i)};
            found = true;
          }
        });
      if (found) { ++update; }
    }, Kokkos::Sum<unsigned>(numValidEntityKeys)
  );
  Kokkos::fence();

  Kokkos::View<Entity*> compactedEntities("locallyAvailEntities", numValidEntityKeys);
  Kokkos::View<EntityWrapper<EntityKey>*> compactedEntityKeys("locallyAvailEntityKeys", numValidEntityKeys);
  Kokkos::parallel_scan(numEntityKeys,
    KOKKOS_LAMBDA(const int j, unsigned& update, bool isFinal) {
      bool isValid = tempEntities(j).is_local_offset_valid();
      unsigned pos = update;

      if (isValid) { ++update; }
      if (isFinal && isValid) {
        compactedEntities(pos) = tempEntities(j);
        compactedEntityKeys(pos) = wrappedEntityKeys(j);
      }
    }
  );
  Kokkos::fence();

  return LocallyAvailEntityKeys{compactedEntityKeys, compactedEntities};
}

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

template<typename ViewType>
ViewType get_sorted_view(const ViewType& view)
{
  using ExecSpace = typename ViewType::execution_space;
  if (!Kokkos::Experimental::is_sorted(ExecSpace{}, view)) {
    ViewType copy("copy_view", view.size());
    Kokkos::deep_copy(ExecSpace{},copy, view);
    Kokkos::sort(ExecSpace{},copy);
    return copy;
  }

  return view;
}

template<class InputIt1, class OutputIt>
KOKKOS_INLINE_FUNCTION
OutputIt my_copy(InputIt1 first, InputIt1 last, OutputIt dest)
{
  for(; first != last; ++first) {
    *dest++ = *first;
  }
  return dest;
}

template<class InputIt1, class InputIt2, class OutputIt>
KOKKOS_INLINE_FUNCTION
OutputIt my_merge(InputIt1 first1, InputIt1 last1,
                  InputIt2 first2, InputIt2 last2,
                  OutputIt d_first)
{
  for (; first1 != last1; ++d_first) {
    if (first2 == last2) {
      return my_copy(first1, last1, d_first);
    }

    if (*first2 < *first1) {
      *d_first = *first2;
      ++first2;
    }
    else {
      *d_first = *first1;
      ++first1;
    }
  }
  return my_copy(first2, last2, d_first);
}

template<typename MESH_TYPE, typename EntityViewType, typename AddPartsViewType, typename NewPartsViewType, typename PartOrdinalsProxyViewType>
void set_add_part_list_per_entity(const MESH_TYPE&,
                                  stk::topology::rank_t rank,
                                  const EntityViewType& entities,
                                  const AddPartsViewType& addParts,
                                  NewPartsViewType& newPartsPerEntity,
                                  PartOrdinalsProxyViewType& partOrdinalsProxyView)
{
  using ExecSpace = typename MESH_TYPE::MeshExecSpace;
  using TeamMember = typename stk::ngp::TeamPolicy<ExecSpace>::member_type;

  STK_NGP_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, addParts));
  auto teamPolicy = stk::ngp::TeamPolicy<ExecSpace>(entities.size(), Kokkos::AUTO);

  Kokkos::parallel_for("set_new_part_lists", teamPolicy,
    KOKKOS_LAMBDA(const TeamMember& teamMember) {
      const unsigned i = teamMember.league_rank();
      unsigned myNumParts = addParts.size();
      const unsigned myStartIdx = i * myNumParts;

      Ordinal* dest = &newPartsPerEntity(myStartIdx);
      Kokkos::single(Kokkos::PerTeam(teamMember),[&]() {
        for (unsigned j = 0; j < myNumParts; ++j) {
          newPartsPerEntity(myStartIdx + j) = addParts(j);
        }
        typename PartOrdinalsProxyViewType::value_type proxyIndices(rank, dest, myNumParts);
        partOrdinalsProxyView(i) = proxyIndices;
      });

      auto begin = Kokkos::Experimental::begin(newPartsPerEntity) + myStartIdx;
      auto end = begin + myNumParts;

      Kokkos::single(Kokkos::PerTeam(teamMember),[&]() {
        partOrdinalsProxyView(i).length = static_cast<unsigned>(Kokkos::Experimental::distance(begin, end));
      });
    }
  );

  Kokkos::fence();
}

template<typename MESH_TYPE, typename EntityViewType, typename AddPartsViewType, typename RmPartsViewType,
         typename NewPartsViewType, typename PartOrdinalsProxyViewType>
void set_new_part_list_per_entity(const MESH_TYPE& ngpMesh,
                                  const EntityViewType& entities,
                                  const AddPartsViewType& addParts,
                                  const RmPartsViewType& rmParts,
                                  unsigned maxPartsPerEntity,
                                  NewPartsViewType& newPartsPerEntity,
                                  PartOrdinalsProxyViewType& partOrdinalsProxyView)
{
  using BucketType = typename MESH_TYPE::BucketType;
  using ExecSpace = typename MESH_TYPE::MeshExecSpace;
  using TeamMember = typename stk::ngp::TeamPolicy<ExecSpace>::member_type;

  STK_NGP_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, addParts));
  auto teamPolicy = stk::ngp::TeamPolicy<ExecSpace>(entities.size(), Kokkos::AUTO);

  Kokkos::parallel_for("set_new_part_lists", teamPolicy,
    KOKKOS_LAMBDA(const TeamMember& teamMember) {
      const unsigned i = teamMember.league_rank();
      const unsigned myStartIdx = i*maxPartsPerEntity;
      const stk::mesh::EntityRank rank = ngpMesh.entity_rank(entities(i));
      const stk::mesh::FastMeshIndex entityIdx = ngpMesh.device_mesh_index(entities(i));
      const BucketType& bucket = ngpMesh.get_bucket(rank, entityIdx.bucket_id);
      auto currentPartOrdsPair = bucket.superset_part_ordinals();

      Ordinal* dest = &newPartsPerEntity(myStartIdx);
      const Ordinal* first1 = addParts.data();
      const Ordinal* last1 = first1+addParts.size();
      const Ordinal* first2 = currentPartOrdsPair.first;
      const Ordinal* last2 = currentPartOrdsPair.second;
      unsigned myNumParts = (last2-first2) + addParts.size();
      Kokkos::single(Kokkos::PerTeam(teamMember),[&]() {
        my_merge(first1, last1, first2, last2, dest);

        typename PartOrdinalsProxyViewType::value_type proxyIndices(rank, dest, myNumParts);
        partOrdinalsProxyView(i) = proxyIndices;
      });

      auto isInRmParts = [&](Ordinal item) {
                           for(unsigned rp=0; rp<rmParts.size(); ++rp)  {
                             if (item == rmParts(rp)) { return true; }
                           }
                           return false;
                         };

      auto begin = Kokkos::Experimental::begin(newPartsPerEntity) + myStartIdx;
      auto end = begin + myNumParts;

      end = Kokkos::Experimental::remove_if(teamMember, begin, end, isInRmParts);
      end = Kokkos::Experimental::unique(teamMember, begin, end);

      Kokkos::single(Kokkos::PerTeam(teamMember),[&]() {
        partOrdinalsProxyView(i).length = static_cast<unsigned>(Kokkos::Experimental::distance(begin, end));
      });
    }
  );

  Kokkos::fence();
}

template <typename ViewType, typename ExecSpace>
void sort_and_unique(ViewType& view, ExecSpace const& execSpace)
{
  Kokkos::sort(view);
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, view));
  Kokkos::Experimental::unique(execSpace, view);
}

template <typename ViewType, typename ExecSpace>
void unique_and_resize(ViewType& view, ExecSpace const& execSpace)
{
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, view));

  auto newEnd = Kokkos::Experimental::unique(execSpace, view);
  auto begin = Kokkos::Experimental::begin(view);
  size_t newSize = Kokkos::Experimental::distance(begin, newEnd);
  if (newSize != view.extent(0)) {
    Kokkos::resize(view, newSize);
  }
}

template <typename ViewType, typename ExecSpace>
void sort_and_unique_and_resize(ViewType& view, ExecSpace const& execSpace)
{
  Kokkos::sort(view);
  unique_and_resize(view, execSpace);
}

template <typename DeviceBucketRepoType, typename PartOrdinalsViewType>
bool has_ranked_part(DeviceBucketRepoType const& deviceBucketRepo, PartOrdinalsViewType const& addPartOrdinals, PartOrdinalsViewType const& removePartOrdinals)
{
  int numRankedParts = 0;

  Kokkos::parallel_reduce(addPartOrdinals.extent(0),
    KOKKOS_LAMBDA(const int& i, int& update) {
      if (deviceBucketRepo.does_induce(addPartOrdinals(i))) {
        update++;
      }
    }, Kokkos::Sum<int>(numRankedParts)
  );
  Kokkos::fence();

  if (numRankedParts > 0) { return true; }

  Kokkos::parallel_reduce(removePartOrdinals.extent(0),
    KOKKOS_LAMBDA(const int& i, int& update) {
      if (deviceBucketRepo.does_induce(removePartOrdinals(i))) {
        update++;
      }
    }, Kokkos::Sum<int>(numRankedParts)
  );
  Kokkos::fence();

  return (numRankedParts > 0);
}

template <typename DeviceBucketRepoType, typename PartOrdinalsViewType>
bool has_ranked_remove_part(DeviceBucketRepoType const& deviceBucketRepo, PartOrdinalsViewType const& removePartOrdinals)
{
  int numRankedParts = 0;

  Kokkos::parallel_reduce(removePartOrdinals.extent(0),
    KOKKOS_LAMBDA(const int& i, int& update) {
      if (deviceBucketRepo.is_ranked_part(removePartOrdinals(i))) {
        update++;
      }
    }, Kokkos::Sum<int>(numRankedParts)
  );
  Kokkos::fence();

  return (numRankedParts > 0);
}

template <typename MeshType, typename EntityViewType>
int get_max_num_downward_connected_entities(MeshType const& ngpMesh, EntityViewType const& entities)
{
  int max = 0;
  auto numEntities = entities.extent(0);

  Kokkos::parallel_reduce(numEntities,
    KOKKOS_LAMBDA (const int& i, int& update) {
      EntityRank rank = ngpMesh.entity_rank(entities(i));

      if (rank == stk::topology::NODE_RANK) {
        update = (update > 0) ? update : 0;
        return;
      }

      auto entityIdx = ngpMesh.device_mesh_index(entities(i));
      auto& bucket = ngpMesh.get_bucket(rank, entityIdx.bucket_id);
      auto totalNumConnectedEntities = 0;

      for (auto rankToSearch = stk::topology::NODE_RANK; rankToSearch < rank; ++rankToSearch) {
        totalNumConnectedEntities += bucket.get_connected_entities(entityIdx.bucket_ord, rankToSearch).size();
      }

      update = (update > totalNumConnectedEntities) ? update : totalNumConnectedEntities;
    }, Kokkos::Max<int>(max)
  );
  Kokkos::fence();

  return max;
}

template <typename WrapperUViewType, typename ExecSpace>
void merge_sorted_duplicate_wrapper_flags(WrapperUViewType const& uview, size_t length, ExecSpace const& execSpace)
{
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(execSpace, 0, length),
    KOKKOS_LAMBDA(const int i) {
      const bool isRunLeader = (i == 0) || !(uview(i).value == uview(i-1).value);
      if (!isRunLeader) { return; }

      bool mergedUserInput = false;
      bool mergedInduction = false;
      int k = i;
      for (; k < static_cast<int>(length) && uview(k).value == uview(i).value; ++k) {
        mergedUserInput = mergedUserInput || uview(k).isUserInputEntity;
        mergedInduction = mergedInduction || uview(k).isForPartInduction;
      }
      for (int j = i; j < k; ++j) {
        uview(j).isUserInputEntity = mergedUserInput;
        uview(j).isForPartInduction = mergedInduction;
      }
    }
  );
  Kokkos::fence();
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

template <typename MeshType, typename EntityViewType, typename EntityKeyViewType>
void populate_applied_entity_keys_for_inducible_part_change(MeshType const& ngpMesh, EntityViewType entities,
                                                             int entityInterval, EntityKeyViewType entityKeys,
                                                             bool hasInduciblePart)
{
  using EntityKeyWrapper = EntityWrapper<EntityKey>;

  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int i) {
      auto myStartIdx = i * entityInterval;
      auto myCurrentIdx = myStartIdx;
      auto entity = entities(i);
      auto rank = ngpMesh.entity_rank(entity);
      auto fastMeshIdx = ngpMesh.device_mesh_index(entity);
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIdx.bucket_id);

      if (!bucket.is_owned()) { return; }

      auto key = ngpMesh.entity_key(entity);
      entityKeys(myCurrentIdx++) = EntityKeyWrapper(key, true);

      if (hasInduciblePart) {
        for (auto iRank = rank; iRank >= stk::topology::NODE_RANK; --iRank) {
          if (iRank == rank) { continue; }
          auto connectedEntities = bucket.get_connected_entities(fastMeshIdx.bucket_ord, iRank);

          for (unsigned j = 0; j < connectedEntities.size(); j++) {
            auto connKey = ngpMesh.entity_key(connectedEntities[j]);
            entityKeys(myCurrentIdx++) = EntityKeyWrapper(connKey, false, true);
          }
        }
      }
    }
  );
  Kokkos::fence();
}

template <typename EntityKeyWrapperViewType, typename MeshType, typename EntityViewType>
EntityKeyWrapperViewType populate_applied_entity_keys(MeshType const& ngpMesh, EntityViewType entities, bool hasInduciblePart)
{
  auto maxNumDownwardConnectedEntities = impl::get_max_num_downward_connected_entities(ngpMesh, entities);
  auto entityInterval = maxNumDownwardConnectedEntities + 1;
  auto maxNumAppliedEntities = entityInterval * entities.extent(0);

  EntityKeyWrapperViewType allAppliedKeys("allAppliedKeys", maxNumAppliedEntities);
  populate_applied_entity_keys_for_inducible_part_change(ngpMesh, entities, entityInterval, allAppliedKeys, hasInduciblePart);

  impl::remove_invalid_entity_keys_sort_unique_and_resize(allAppliedKeys, typename MeshType::MeshExecSpace{});

  return allAppliedKeys;
}

template <typename MeshType, typename EntityKeyWrapperViewType>
EntityKeyWrapperViewType intersect_applied_entity_keys_with_comm_map(MeshType const& ngpMesh,
                                                                     EntityKeyWrapperViewType const& appliedEntityKeys,
                                                                     int destProc, bool includeGhosts = true)
{
  using ExecSpace = typename MeshType::MeshExecSpace;

  auto numAppliedEntityKeys = appliedEntityKeys.extent(0);
  if (numAppliedEntityKeys == 0) {
    return EntityKeyWrapperViewType("intersectedEntityKeys", 0);
  }

  Kokkos::View<int*> keysInIntersection("keysInIntersection", numAppliedEntityKeys);

  unsigned numEntityKeysFound = 0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, numAppliedEntityKeys),
    KOKKOS_LAMBDA(const int i, unsigned& update) {
      auto key = appliedEntityKeys(i).value;
      auto rank = key.rank();

      bool isRelevantToDestProc = find_entity_key_in_comm_map(ngpMesh, rank, destProc, key,
                                                              includeGhosts).is_local_offset_valid();
      if (isRelevantToDestProc) {
        keysInIntersection(i) = 1;
        ++update;
      }
    }, Kokkos::Sum<unsigned>(numEntityKeysFound)
  );
  Kokkos::fence();

  EntityKeyWrapperViewType intersectedKeys("intersectedKeys", numEntityKeysFound);
  if (numEntityKeysFound == 0) { return intersectedKeys; }

  Kokkos::parallel_scan(Kokkos::RangePolicy<ExecSpace>(0, numAppliedEntityKeys),
    KOKKOS_LAMBDA(const size_t i, unsigned& update, const bool isFinal) {
      bool isFoundKey = (keysInIntersection(i) != 0);
      unsigned pos = update;

      if (isFoundKey) { ++update; }
      if (isFinal && isFoundKey) { intersectedKeys(pos) = appliedEntityKeys(i); }
    }
  );
  Kokkos::fence();

  return intersectedKeys;
}

template <typename MeshType, typename InputWrappedEntityViewType, typename WrappedEntityViewType>
void populate_all_downward_connected_entities_and_wrap_entities(MeshType const& ngpMesh, InputWrappedEntityViewType const& entities,
                                                                int entityInterval, WrappedEntityViewType const& wrappedEntities)
{
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int i) {
      auto myStartIdx = i * entityInterval;
      auto myCurrentIdx = myStartIdx;
      auto rank = ngpMesh.entity_rank(entities(i));
      auto fastMeshIdx = ngpMesh.device_mesh_index(entities(i));
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIdx.bucket_id);

      auto wrappedEntity = entities(i);
      if (wrappedEntity.isUserInputEntity || wrappedEntity.isForPartInduction)
        wrappedEntities(myCurrentIdx++) = EntityWrapper(wrappedEntity.value, wrappedEntity.isUserInputEntity, wrappedEntity.isForPartInduction);
      else
        wrappedEntities(myCurrentIdx++) = EntityWrapper(wrappedEntity.value, true);

      for (auto iRank = rank; iRank >= stk::topology::NODE_RANK; --iRank) {
        if (iRank == rank) { continue; }
        auto connectedEntities = bucket.get_connected_entities(fastMeshIdx.bucket_ord, iRank);

        for (unsigned j = 0; j < connectedEntities.size(); j++) {
          wrappedEntity = EntityWrapper(connectedEntities[j], false, true);
          wrappedEntities(myCurrentIdx++) = wrappedEntity;
        }
      }
    }
  );
  Kokkos::fence();
}

template <typename EntityViewType, typename ExecSpace>
void remove_invalid_entities_sort_unique_and_resize(EntityViewType& entities, ExecSpace const& execSpace)
{
  using EntityType = typename EntityViewType::value_type;
  using MemorySpace = typename EntityViewType::memory_space;
  using EntityUViewType = Kokkos::View<EntityType*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  auto isInvalidEntity = KOKKOS_LAMBDA(EntityType entity) {
                              return static_cast<stk::mesh::Entity>(entity).m_value == 0;
                         };
  auto end = Kokkos::Experimental::remove_if(execSpace, entities, isInvalidEntity);

  auto begin = Kokkos::Experimental::begin(entities);
  auto length = Kokkos::Experimental::distance(begin, end);

  EntityUViewType uview(entities.data(), length);
  Kokkos::sort(uview);
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, uview));

  auto newEnd = Kokkos::Experimental::unique(execSpace, uview);
  auto newBegin = Kokkos::Experimental::begin(uview);
  size_t newSize = Kokkos::Experimental::distance(newBegin, newEnd);
  if (newSize != entities.extent(0)) {
    Kokkos::resize(entities, newSize);
  }
}

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

template <typename EntityViewType, typename ExecSpace>
void remove_invalid_wrapped_entities_sort_unique_merge_and_resize(EntityViewType& wrappedEntities, ExecSpace const& execSpace)
{
  using EntityType = typename EntityViewType::value_type;
  using MemorySpace = typename EntityViewType::memory_space;
  using EntityUViewType = Kokkos::View<EntityType*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  auto isInvalidEntity = KOKKOS_LAMBDA(EntityType entity) {
                              return static_cast<stk::mesh::Entity>(entity).m_value == 0;
                         };
  auto end = Kokkos::Experimental::remove_if(execSpace, wrappedEntities, isInvalidEntity);

  auto begin = Kokkos::Experimental::begin(wrappedEntities);
  auto length = Kokkos::Experimental::distance(begin, end);

  EntityUViewType uview(wrappedEntities.data(), length);
  Kokkos::sort(uview);
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, uview));

  merge_sorted_duplicate_wrapper_flags(uview, static_cast<size_t>(length), execSpace);

  auto newEnd = Kokkos::Experimental::unique(execSpace, uview);
  auto newBegin = Kokkos::Experimental::begin(uview);
  size_t newSize = Kokkos::Experimental::distance(newBegin, newEnd);
  if (newSize != wrappedEntities.extent(0)) {
    Kokkos::resize(wrappedEntities, newSize);
  }
}

template <typename MeshType, typename TeamMember, typename BucketType, typename IsLosingPartFn>
KOKKOS_INLINE_FUNCTION
bool should_keep_induced_part(MeshType const& ngpMesh, TeamMember const& team,
                              BucketType const& bucket, unsigned bucketOrd,
                              PartOrdinal partOrdinal, stk::mesh::EntityRank partRank,
                              IsLosingPartFn const& isLosingPart, bool ownedInducersOnly)
{
  auto connectedUpperRankEntities = bucket.get_connected_entities(bucketOrd, partRank);
  bool anyKeeps = false;

  Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, connectedUpperRankEntities.size()),
    [&](const int& connIdx, bool& localKeeps) {
      auto upperEntity = connectedUpperRankEntities[connIdx];
      auto bucketId = ngpMesh.device_mesh_index(upperEntity).bucket_id;
      auto& upperBucket = ngpMesh.get_bucket(partRank, bucketId);

      if (ownedInducersOnly && !upperBucket.is_owned()) { return; }
      if (!upperBucket.member(partOrdinal)) { return; }
      if (!isLosingPart(upperEntity)) { localKeeps = true; }
    }, Kokkos::LOr<bool>(anyKeeps)
  );

  return anyKeeps;
}

template <typename DeniedPartRemovalViewType>
KOKKOS_INLINE_FUNCTION
bool is_part_removal_denied(DeniedPartRemovalViewType deniedPartRemovalView, Entity entity, PartOrdinal part)
{
  const size_t numDenied = deniedPartRemovalView.extent(0);
  if (numDenied == 0) { return false; }

  const DeniedPartRemoval target{entity, part};
  size_t lo = 0;
  size_t hi = numDenied;
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    if (deniedPartRemovalView(mid) < target) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return (lo < numDenied) && (deniedPartRemovalView(lo) == target);
}

template<typename MeshType, typename WrappedEntityViewType, typename AddPartsViewType, typename RmPartsViewType,
         typename NewPartsViewType, typename PartOrdinalsProxyViewType,
         typename DeniedPartRemovalViewType>
void set_new_part_list_per_entity_with_induced_parts(MeshType const& ngpMesh,
                                                     WrappedEntityViewType wrappedEntities,
                                                     AddPartsViewType addParts,
                                                     RmPartsViewType rmParts,
                                                     unsigned maxPartsPerEntity,
                                                     NewPartsViewType newPartsPerEntity,
                                                     PartOrdinalsProxyViewType partOrdinalsProxyView,
                                                     DeniedPartRemovalViewType deniedPartRemovalView)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using TeamMember = typename stk::ngp::TeamPolicy<ExecSpace>::member_type;

  STK_NGP_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, addParts));
  auto teamPolicy = stk::ngp::TeamPolicy<ExecSpace>(wrappedEntities.size(), Kokkos::AUTO);
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for("set_new_part_lists_with_induced_parts", teamPolicy,
    KOKKOS_LAMBDA(const TeamMember& team) {
      auto idx = team.league_rank();
      auto myStartIdx = idx * maxPartsPerEntity;
      auto rank = ngpMesh.entity_rank(wrappedEntities(idx));
      auto fastMeshIdx = ngpMesh.device_mesh_index(wrappedEntities(idx));
      auto isForPartInduction = wrappedEntities(idx).isForPartInduction;
      auto isUserInputEntity = wrappedEntities(idx).isUserInputEntity;
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIdx.bucket_id);
      auto currentPartOrdsPair = bucket.superset_part_ordinals();

      Ordinal* dest = &newPartsPerEntity(myStartIdx);
      const Ordinal* first1 = addParts.data();
      const Ordinal* last1 = first1+addParts.size();
      const Ordinal* first2 = currentPartOrdsPair.first;
      const Ordinal* last2 = currentPartOrdsPair.second;
      unsigned myNumParts = (last2-first2) + addParts.size();

      Kokkos::single(Kokkos::PerTeam(team),[&]() {
        my_merge(first1, last1, first2, last2, dest);

        typename PartOrdinalsProxyViewType::value_type proxyIndices(rank, dest, myNumParts);
        partOrdinalsProxyView(idx) = proxyIndices;
      });
      team.team_barrier();

      team.team_barrier();

      auto isRemovedPart = [&](Ordinal partOrdinal) { return partOrdinal == InvalidPartOrdinal; };

      auto isInducingPart = [&](Ordinal partOrdinal) { return deviceBucketRepo.does_induce(partOrdinal); };

      auto isInAddParts =  [&](Ordinal partOrdinal) {
                             for (unsigned ap = 0; ap < addParts.size(); ++ap)  {
                               if (partOrdinal == addParts(ap)) { return true; }
                             }
                             return false;
                           };

      auto isInRmParts =  [&](Ordinal partOrdinal) {
                            for (unsigned rp = 0; rp < rmParts.size(); ++rp)  {
                              if (partOrdinal == rmParts(rp)) { return true; }
                            }
                            return false;
                          };
      
      auto isInWrappedEntities = [&](Entity entity) {
                                   for (unsigned re = 0; re < wrappedEntities.extent(0); ++re) {
                                     if (entity == static_cast<Entity>(wrappedEntities(re))) {
                                       return true;
                                     }
                                   }
                                   return false;
                                 };

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, myNumParts),
        [&](const int& thIdx) {
          auto myPartOrdinalIdx = myStartIdx + thIdx;
          auto partToCheck = newPartsPerEntity(myPartOrdinalIdx);
          auto partToCheckRank = deviceBucketRepo.get_part_rank(partToCheck);

          auto inAddParts = isInAddParts(partToCheck);
          auto inRmParts = isInRmParts(partToCheck);
          auto inducingPart = isInducingPart(partToCheck);

          if (!inducingPart && inAddParts) {
            if (!isUserInputEntity) {
              newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
              return;
            }
          } else if(!inducingPart && inRmParts) {
            if (isUserInputEntity) {
              newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
              return;
            }
          }
          else if (inducingPart && inAddParts) {
            if (rank == partToCheckRank) {
              return;
            } else if (partToCheckRank > rank) {
              if (!isForPartInduction) {
                newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
                return;
              }
            }
          } else if (inducingPart && inRmParts) {
            if (rank == partToCheckRank) {
              newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
              return;
            } else if (partToCheckRank > rank) {
              auto isLosingPart = [&](Entity upperEntity) { return isInWrappedEntities(upperEntity); };
              bool keepsPart = should_keep_induced_part(ngpMesh, team, bucket, 0,
                                                        partToCheck, partToCheckRank,
                                                        isLosingPart, false);

              Entity downwardEntity = static_cast<Entity>(wrappedEntities(idx));
              Kokkos::single(Kokkos::PerThread(team), [&]() {
                if (!keepsPart && !is_part_removal_denied(deniedPartRemovalView, downwardEntity, partToCheck)) {
                  newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
                }
              });
              return;
            }
          }
        }
      );
      team.team_barrier();

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
          auto inducingPart = deviceBucketRepo.does_induce(partToCheck);

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
          if (deviceBucketRepo.does_induce(partOrdinal) &&
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

template <typename MeshType, typename EntitySrcDestView, typename NewBucketsToAddViewType>
void assign_initial_bucket_id_and_ordinal(MeshType const& ngpMesh, stk::topology::rank_t rank, EntitySrcDestView& entitySrcDestView,
                                          NewBucketsToAddViewType& numNewBucketsToAddInPartitions)
{
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  deviceBucketRepo.check_mesh_consistency();

  auto begin = Kokkos::Experimental::begin(entitySrcDestView);
  auto end = Kokkos::Experimental::end(entitySrcDestView);

  if (!Kokkos::Experimental::is_sorted(typename MeshType::MeshExecSpace{}, entitySrcDestView)) {
    Kokkos::sort(entitySrcDestView);
  }

  // TODO convert to block partitioned parallel scan
  Kokkos::parallel_for(
      1, KOKKOS_LAMBDA(const int) {
        auto currRank = stk::topology::INVALID_RANK;
        unsigned currDestPartitionId = INVALID_PARTITION_ID;
        unsigned nextBucketOrdinalInPartition = INVALID_INDEX;
        unsigned nextEntityOrdinalInBucket = INVALID_INDEX;

        int newBucketPerPartitionIdx = -1;

        unsigned bucketCapacity = deviceBucketRepo.get_bucket_capacity();

        for (unsigned i = 0; i < entitySrcDestView.extent(0); ++i) {
          auto& srcDest = entitySrcDestView(i);
          auto destPartitionId = srcDest.destPartitionId;

          // new sets of entities into another partition
          if (currDestPartitionId != destPartitionId || currRank != rank) {
            currDestPartitionId = destPartitionId;
            currRank = rank;

            auto partition = deviceBucketRepo.get_partition(rank, destPartitionId);
            auto lastBucketIdx = partition->get_last_avail_bucket_index();

            numNewBucketsToAddInPartitions(++newBucketPerPartitionIdx) = {rank, destPartitionId, 0};

            // append to the end of existing bucket
            if (lastBucketIdx != INVALID_INDEX) {
              auto destBucket = partition->m_buckets[lastBucketIdx];
              nextBucketOrdinalInPartition = lastBucketIdx;
              nextEntityOrdinalInBucket = destBucket->size();
            }
            // place into a new bucket
            else {
              nextBucketOrdinalInPartition = partition->m_buckets.size();
              nextEntityOrdinalInBucket = 0;
              numNewBucketsToAddInPartitions(newBucketPerPartitionIdx).numBucketsToAdd++;
            }
          }

          // need to place into next (new) bucket
          if (nextEntityOrdinalInBucket >= bucketCapacity) {
            nextBucketOrdinalInPartition++;
            nextEntityOrdinalInBucket = 0;
            numNewBucketsToAddInPartitions(newBucketPerPartitionIdx).numBucketsToAdd++;
          }

          srcDest.destBucketId = INVALID_BUCKET_ID;  // to be filled in later; can't be determined now
          srcDest.destBucketIndexInPartition = nextBucketOrdinalInPartition;
          srcDest.destBucketOrd = nextEntityOrdinalInBucket;

          nextEntityOrdinalInBucket++;
        }
      });
  Kokkos::fence();
}

template <typename MeshType, typename EntitySrcDestView, typename NewBucketsToAddViewType,
          typename GrowLastBucketViewType>
void assign_dest_bucket_id_and_ordinal(MeshType const& ngpMesh, EntitySrcDestView& entitySrcDestView,
                                       NewBucketsToAddViewType& numNewBucketsToAddInPartitions,
                                       GrowLastBucketViewType& growLastBucketInPartition)
{
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  deviceBucketRepo.check_mesh_consistency();

  auto begin = Kokkos::Experimental::begin(entitySrcDestView);
  auto end = Kokkos::Experimental::remove_if(
      typename MeshType::MeshExecSpace{}, entitySrcDestView,
      KOKKOS_LAMBDA(typename EntitySrcDestView::value_type const& srcDest) {
        return srcDest.srcPartitionId == srcDest.destPartitionId;
      });

  auto newLength = Kokkos::Experimental::distance(begin, end);
  Kokkos::resize(entitySrcDestView, newLength); 

  if (!Kokkos::Experimental::is_sorted(typename MeshType::MeshExecSpace{}, entitySrcDestView)) {
    Kokkos::sort(entitySrcDestView);
  }

  // TODO convert to block partitioned parallel scan
  Kokkos::parallel_for(
      1, KOKKOS_LAMBDA(const int) {
        auto currRank = stk::topology::INVALID_RANK;
        unsigned currDestPartitionId = INVALID_PARTITION_ID;
        unsigned nextBucketOrdinalInPartition = INVALID_INDEX;
        unsigned nextEntityOrdinalInBucket = INVALID_INDEX;

        int newBucketPerPartitionIdx = -1;
        int growLastBucketInPartitionIdx = -1;

        unsigned maxBucketCapacity = deviceBucketRepo.get_maximum_bucket_capacity();

        for (unsigned i = 0; i < entitySrcDestView.extent(0); ++i) {
          auto& srcDest = entitySrcDestView(i);
          auto rank = srcDest.rank;
          auto destPartitionId = srcDest.destPartitionId;
          unsigned destBucketCapacity = maxBucketCapacity;

          // new sets of entities into another partition
          if (currDestPartitionId != destPartitionId || currRank != rank) {
            currDestPartitionId = destPartitionId;
            currRank = rank;

            auto partition = deviceBucketRepo.get_partition(rank, destPartitionId);
            auto lastBucketIdx = partition->get_last_avail_bucket_index(maxBucketCapacity);

            numNewBucketsToAddInPartitions(++newBucketPerPartitionIdx) = {rank, destPartitionId, 0};

            if (lastBucketIdx != INVALID_INDEX) {
              // append to the end of existing bucket
              auto destBucket = partition->m_buckets[lastBucketIdx];
              nextBucketOrdinalInPartition = lastBucketIdx;
              nextEntityOrdinalInBucket = destBucket->size();
              destBucketCapacity = destBucket->capacity();
            }
            else {
              // place into a new bucket
              nextBucketOrdinalInPartition = partition->m_buckets.size();
              nextEntityOrdinalInBucket = 0;
              numNewBucketsToAddInPartitions(newBucketPerPartitionIdx).numBucketsToAdd++;
            }
          }

          if (nextEntityOrdinalInBucket >= maxBucketCapacity) {
            // need to place into next (new) bucket
            nextBucketOrdinalInPartition++;
            nextEntityOrdinalInBucket = 0;
            numNewBucketsToAddInPartitions(newBucketPerPartitionIdx).numBucketsToAdd++;
          }
          else if (nextEntityOrdinalInBucket >= destBucketCapacity) {
            growLastBucketInPartition(++growLastBucketInPartitionIdx) = {rank, destPartitionId, true};
          }

          srcDest.destBucketId = INVALID_BUCKET_ID;  // to be filled in later; can't be determined now
          srcDest.destBucketIndexInPartition = nextBucketOrdinalInPartition;
          srcDest.destBucketOrd = nextEntityOrdinalInBucket;

          // printf("assign_dest_bucket_id_and_ordinal at [%u] entitySrcDest:\n"
          //        "\tsrcPartitionId = %u\n",
          //        i,
          //        srcDest.srcPartitionId);

          // printf("\tdestPartitionId = %u\n"
          //        "\tdestBucketIndexInPartition = %u\n"
          //        "\tdestBucketOrd = %u\n"
          //        "\tnum new buckets to add = %u\n"
          //        "\tdestNodeConnectivityStartIdx = %u\n"
          //        "\tdestNodeConnectivityLength = %u\n"
          //        srcDest.destPartitionId, srcDest.destBucketIndexInPartition, srcDest.destBucketOrd,
          //        numNewBucketsToAddInPartitions(newBucketPerPartitionIdx).numBucketsToAdd,

          nextEntityOrdinalInBucket++;
        }
      });
  Kokkos::fence();
}

template <typename MeshType, typename EntitySrcDestView>
void set_dest_bucket_ids(MeshType const& ngpMesh, EntitySrcDestView& entitySrcDestView)
{
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for(entitySrcDestView.extent(0),
    KOKKOS_LAMBDA(const int i) {
      auto& entitySrcDest = entitySrcDestView(i);
      auto rank = entitySrcDest.rank;
      auto destPartitionId = entitySrcDest.destPartitionId;
      auto destPartition = deviceBucketRepo.get_partition(rank, destPartitionId);
      auto destBucketId = destPartition->get_bucket(entitySrcDest.destBucketIndexInPartition)->bucket_id();
      entitySrcDest.destBucketId = destBucketId;
    }
  );
}

template <typename ExecSpace, typename ViewType, typename T>
ViewType remove_invalid_entries_and_resize(ViewType const& view, T const invalidValue)
{
  auto begin = Kokkos::Experimental::begin(view);
  auto end = Kokkos::Experimental::end(view);

  using value_type = typename ViewType::non_const_value_type;
  auto pred = KOKKOS_LAMBDA(const value_type v) -> bool { return v != invalidValue; };

  auto count = Kokkos::Experimental::count_if(ExecSpace{}, begin, end, pred);

  ViewType filteredView("filteredView", count);
  auto insertBegin = Kokkos::Experimental::begin(filteredView);

  Kokkos::Experimental::copy_if(ExecSpace{}, begin, end, insertBegin, pred);
  return filteredView;
}

template <typename MeshType, typename EntityView>
void require_entity_owner(MeshType& ngpMesh, EntityView entities)
{
  bool anyNonOwned = false;
  Kokkos::parallel_reduce(entities.extent(0),
    KOKKOS_LAMBDA(const int i, bool& update) {
      auto entity = entities(i);
      auto rank = ngpMesh.entity_rank(entity);
      auto fmi = ngpMesh.fast_mesh_index(entity);
      auto& bucket = ngpMesh.get_bucket(rank, fmi.bucket_id);
      auto isOwned = bucket.is_owned();
      update |= !isOwned;
    }, Kokkos::LOr<bool>(anyNonOwned)
  );

  STK_ThrowRequireMsg(!anyNonOwned, "All entities must be owned");
}

template <typename MeshType, typename EntityView, typename PartOrdinalView>
void communicate_shared_entities(MeshType& ngpMesh, CommSparse& commSparse, EntityView entities,
                                 PartOrdinalView addPartOrdinals, PartOrdinalView removePartOrdinals)
{
  using EntityKeyWrapper = EntityWrapper<EntityKey>;
  using EntityKeyWrapperView = Kokkos::View<EntityKeyWrapper*>;
  using EntityKeyWrapperHostView = Kokkos::View<EntityKeyWrapper*, stk::ngp::HostMemSpace>;

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  bool hasRankedPart = has_ranked_part(deviceBucketRepo, addPartOrdinals, removePartOrdinals);
  auto appliedEntityKeys = populate_applied_entity_keys<EntityKeyWrapperView>(ngpMesh, entities, hasRankedPart);
  auto addPartOrdinalsOnHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, addPartOrdinals);
  auto removePartOrdinalsOnHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, removePartOrdinals);

  auto comm = ngpMesh.get_bulk_on_host().parallel();
  auto numProcs = commSparse.parallel_size();
  auto myRank = stk::parallel_machine_rank(comm);

  std::vector<EntityKeyWrapperHostView> entityKeysToSendPerProc(numProcs);
  for (int destProc = 0; destProc < numProcs; ++destProc) {
    if (destProc == myRank) { continue; }

    auto entityKeysForDestProc = intersect_applied_entity_keys_with_comm_map(ngpMesh, appliedEntityKeys, destProc);
    entityKeysToSendPerProc[destProc] = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, entityKeysForDestProc);
  }

  stk::pack_and_communicate(commSparse, [&]() {
    for (int destProc = 0; destProc < numProcs; ++destProc) {
      if (destProc == myRank) { continue; }

      auto& entityKeysToSend = entityKeysToSendPerProc[destProc];
      unsigned numSharedEntityKeys = entityKeysToSend.extent(0);
      unsigned numAddPartOrdinals = addPartOrdinalsOnHost.extent(0);
      unsigned numRemovePartOrdinals = removePartOrdinalsOnHost.extent(0);

      auto& buf = commSparse.send_buffer(destProc);
      buf.pack(numSharedEntityKeys);
      buf.pack(numAddPartOrdinals);
      buf.pack(numRemovePartOrdinals);

      if (numSharedEntityKeys == 0) { continue; }

      for (unsigned i = 0; i < numSharedEntityKeys; ++i)
        buf.pack(entityKeysToSend(i));

      for (unsigned i = 0; i < numAddPartOrdinals; ++i)
        buf.pack(addPartOrdinalsOnHost(i));

      for (unsigned i = 0; i < numRemovePartOrdinals; ++i)
        buf.pack(removePartOrdinalsOnHost(i));
    }
  });
}

template <typename MeshType, typename CallbackFn>
void unpack_shared_entities_and_callback(MeshType& ngpMesh, CommSparse& commSparse, CallbackFn callback)
{
  using EntityKeyWrapper = EntityWrapper<EntityKey>;

  auto comm = stk::parallel_machine_world();
  auto numProcs = commSparse.parallel_size();
  auto myRank = stk::parallel_machine_rank(comm);

  std::vector<EntityKeyWrapper> entityKeys;
  std::vector<unsigned> addPartOrdinals;
  std::vector<unsigned> removePartOrdinals;

  for (int srcProc = 0; srcProc < numProcs; ++srcProc) {
    if (srcProc == myRank) { continue; }

    unsigned numSharedEntityKeys = 0;
    unsigned numAddPartOrdinals = 0;
    unsigned numRemovePartOrdinals = 0;

    stk::CommBuffer& buf = commSparse.recv_buffer(srcProc);
    buf.unpack(numSharedEntityKeys);
    buf.unpack(numAddPartOrdinals);
    buf.unpack(numRemovePartOrdinals);

    if (numSharedEntityKeys == 0) { continue; }

    auto unpack_vals = [&]<typename UnpackType>(auto& unpackCount, auto& unpackTo) {
      unpackTo.clear();
      for (unsigned i = 0; i < unpackCount; ++i) {
        UnpackType val;
        buf.unpack(val);
        unpackTo.push_back(val);
      }
    };

    unpack_vals.template operator()<EntityKeyWrapper>(numSharedEntityKeys, entityKeys);
    unpack_vals.template operator()<unsigned>(numAddPartOrdinals, addPartOrdinals);
    unpack_vals.template operator()<unsigned>(numRemovePartOrdinals, removePartOrdinals);

    using EntityKeysHostView = Kokkos::View<EntityKeyWrapper*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using PartOrdinalsHostView = Kokkos::View<unsigned*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    EntityKeysHostView sharedEntityKeysOnHost(entityKeys.data(), entityKeys.size());
    auto sharedEntityKeysView = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, sharedEntityKeysOnHost);

    auto locallyAvailEntityKeys = get_locally_avail_entity_keys(ngpMesh, sharedEntityKeysView);

    if (locallyAvailEntityKeys.entities.extent(0) == 0) { continue; }

    auto wrappedSharedEntityView = convert_to_wrapped_entities(locallyAvailEntityKeys.entities, locallyAvailEntityKeys.keys);

    PartOrdinalsHostView addPartOrdinalsOnHost(addPartOrdinals.data(), addPartOrdinals.size());
    Kokkos::View<unsigned*> addPartOrdinalsView("", addPartOrdinals.size());
    Kokkos::deep_copy(addPartOrdinalsView, addPartOrdinalsOnHost);

    PartOrdinalsHostView removePartOrdinalsOnHost(removePartOrdinals.data(), removePartOrdinals.size());
    Kokkos::View<unsigned*> removePartOrdinalsView("", removePartOrdinals.size());
    Kokkos::deep_copy(removePartOrdinalsView, removePartOrdinalsOnHost);

    callback(wrappedSharedEntityView, addPartOrdinalsView, removePartOrdinalsView);
  }
}

inline void unpack_communicated_shared_entities(stk::CommBuffer& buf,
                                                std::vector<EntityWrapper<EntityKey>>& entityKeys,
                                                std::vector<unsigned>& removePartOrdinals)
{
  unsigned numSharedEntityKeys = 0;
  unsigned numAddPartOrdinals = 0;
  unsigned numRemovePartOrdinals = 0;
  buf.unpack(numSharedEntityKeys);
  buf.unpack(numAddPartOrdinals);
  buf.unpack(numRemovePartOrdinals);

  if (numSharedEntityKeys != 0) {
    for (unsigned i = 0; i < numSharedEntityKeys; ++i) {
      EntityWrapper<EntityKey> entityKey;
      buf.unpack(entityKey);
      entityKeys.push_back(entityKey);
    }
    for (unsigned i = 0; i < numAddPartOrdinals; ++i) {
      unsigned addOrdinal;
      buf.unpack(addOrdinal);
    }
    for (unsigned i = 0; i < numRemovePartOrdinals; ++i) {
      unsigned removeOrdinal;
      buf.unpack(removeOrdinal);
      removePartOrdinals.push_back(removeOrdinal);
    }
  }

  buf.reset();
}

template <typename MeshType, typename EntityView, typename PartOrdinalView>
Kokkos::View<int*, stk::ngp::HostMemSpace>
compute_denied_removal_flags(MeshType& ngpMesh,
                             std::vector<EntityWrapper<EntityKey>>& remoteSharedEntityKeys,
                             EntityView localInputEntities, PartOrdinalView localInputRemovePartOrdinals,
                             std::vector<unsigned>& remoteRemovePartOrdinals,
                             std::vector<EntityWrapper<EntityKey>>& compactedRemoteEntityKeys,
                             size_t& numSharedEntitiesWithSrcProc)
{
  using EntityKeyWrapper = EntityWrapper<EntityKey>;
  using ExecSpace = typename MeshType::MeshExecSpace;
  using TeamMember = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  using EntityKeysHostView = Kokkos::View<EntityKeyWrapper*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using PartOrdinalsHostView = Kokkos::View<unsigned*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  EntityKeysHostView remoteSharedEntityKeysHost(remoteSharedEntityKeys.data(), remoteSharedEntityKeys.size());
  auto remoteSharedEntityKeysDevice = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, remoteSharedEntityKeysHost);
  auto locallyAvailEntityKeys = get_locally_avail_entity_keys(ngpMesh, remoteSharedEntityKeysDevice);
  auto remoteSharedEntities = locallyAvailEntityKeys.entities;

  numSharedEntitiesWithSrcProc = remoteSharedEntities.extent(0);

  auto compactedKeysHost = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, locallyAvailEntityKeys.keys);
  compactedRemoteEntityKeys.assign(compactedKeysHost.data(), compactedKeysHost.data() + compactedKeysHost.extent(0));

  auto numRemoteRemoveParts = remoteRemovePartOrdinals.size();

  Kokkos::View<int*> deniedFlags("deniedFlags", numSharedEntitiesWithSrcProc * numRemoteRemoveParts);
  if (numSharedEntitiesWithSrcProc == 0) {
    return Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, deniedFlags);
  }

  PartOrdinalsHostView remoteRemovePartOrdinalsHost(remoteRemovePartOrdinals.data(), remoteRemovePartOrdinals.size());
  Kokkos::View<unsigned*> remoteRemovePartOrdinalsDevice("remoteRemovePartOrdinalsDevice", remoteRemovePartOrdinals.size());
  Kokkos::deep_copy(remoteRemovePartOrdinalsDevice, remoteRemovePartOrdinalsHost);

  auto teamPolicy = Kokkos::TeamPolicy<ExecSpace>(numSharedEntitiesWithSrcProc, Kokkos::AUTO);

  Kokkos::parallel_for("compute_denied_part_removals", teamPolicy,
    KOKKOS_LAMBDA(TeamMember const& team) {
      auto entityIndex = team.league_rank();
      Entity downwardEntity = remoteSharedEntities(entityIndex);
      auto rank = ngpMesh.entity_rank(downwardEntity);
      auto fmi = ngpMesh.device_mesh_index(downwardEntity);
      auto& bucket = ngpMesh.get_bucket(rank, fmi.bucket_id);

      for (size_t removeIndex = 0; removeIndex < numRemoteRemoveParts; ++removeIndex) {
        auto partToCheck = remoteRemovePartOrdinalsDevice(removeIndex);
        if (!deviceBucketRepo.is_ranked_part(partToCheck)) { continue; }
        auto partToCheckRank = deviceBucketRepo.get_part_rank(partToCheck);
        if (!(partToCheckRank > rank)) { continue; }

        bool isInInputRemovesPart = false;
        for (unsigned k = 0; k < localInputRemovePartOrdinals.extent(0); ++k) {
          if (localInputRemovePartOrdinals(k) == partToCheck) {
            isInInputRemovesPart = true;
            break;
          }
        }

        auto isLosingPart = [&](Entity upperEntity) {
          if (!isInInputRemovesPart) { return false; }

          for (unsigned k = 0; k < localInputEntities.extent(0); ++k) {
            if (localInputEntities(k) == upperEntity) { return true; }
          }
          return false;
        };

        bool shouldKeepPart = should_keep_induced_part(ngpMesh, team, bucket, fmi.bucket_ord,
                                                       partToCheck, partToCheckRank,
                                                       isLosingPart, true);

        Kokkos::single(Kokkos::PerTeam(team), [&]() {
          if (shouldKeepPart) {
            deniedFlags(entityIndex * numRemoteRemoveParts + removeIndex) = 1;
          }
        });
      }
    }
  );
  Kokkos::fence();

  return Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, deniedFlags);
}

template <typename MeshType, typename EntityView, typename PartOrdinalView>
void compute_denied_parts_to_remove(MeshType& ngpMesh, CommSparse& commSparse,
                                    EntityView localEntities, PartOrdinalView localRemovePartOrdinals,
                                    std::vector<std::vector<EntityWrapper<EntityKey>>>& deniedKeysPerProc,
                                    std::vector<std::vector<unsigned>>& deniedRemovalPartsPerProc)
{
  using EntityKeyWrapper = EntityWrapper<EntityKey>;

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  auto comm = commSparse.parallel();
  auto numProcs = commSparse.parallel_size();
  auto myRank = stk::parallel_machine_rank(comm);

  for (int srcProc = 0; srcProc < numProcs; ++srcProc) {
    if (srcProc == myRank) { continue; }

    std::vector<EntityKeyWrapper> remoteSharedEntityKeys;
    std::vector<unsigned> remotePartOrdsToRemove;
    unpack_communicated_shared_entities(commSparse.recv_buffer(srcProc), remoteSharedEntityKeys, remotePartOrdsToRemove);

    if (remoteSharedEntityKeys.empty() || remotePartOrdsToRemove.empty()) { continue; }

    bool anyRankedRemovePart = false;
    for (unsigned removeOrdinal : remotePartOrdsToRemove) {
      if (deviceBucketRepo.is_ranked_part(removeOrdinal)) {
        anyRankedRemovePart = true;
        break;
      }
    }
    if (!anyRankedRemovePart) { continue; }

    size_t numSharedEntitiesWithSrcProc = 0;
    std::vector<EntityKeyWrapper> compactedRemoteEntityKeys;
    auto deniedFlags = compute_denied_removal_flags(ngpMesh, remoteSharedEntityKeys, localEntities, localRemovePartOrdinals,
                                                   remotePartOrdsToRemove, compactedRemoteEntityKeys,
                                                   numSharedEntitiesWithSrcProc);

    auto numRemoveParts = remotePartOrdsToRemove.size();
    for (size_t entityIndex = 0; entityIndex < numSharedEntitiesWithSrcProc; ++entityIndex) {
      for (size_t removeIndex = 0; removeIndex < numRemoveParts; ++removeIndex) {
        if (deniedFlags(entityIndex * numRemoveParts + removeIndex) != 0) {
          deniedKeysPerProc[srcProc].push_back(compactedRemoteEntityKeys[entityIndex]);
          deniedRemovalPartsPerProc[srcProc].push_back(remotePartOrdsToRemove[removeIndex]);
        }
      }
    }
  }
}

template <typename MeshType, typename EntityView, typename PartOrdinalView>
void communicate_denied_part_removals(MeshType& ngpMesh, CommSparse& commSparse, CommSparse& partRemovalDenialCommSparse,
                                      EntityView localEntities, PartOrdinalView localRemovePartOrdinals)
{
  using EntityKeyWrapper = EntityWrapper<EntityKey>;

  auto comm = commSparse.parallel();
  auto numProcs = commSparse.parallel_size();
  auto myRank = stk::parallel_machine_rank(comm);

  std::vector<std::vector<EntityKeyWrapper>> deniedKeysPerProc(numProcs);
  std::vector<std::vector<unsigned>> deniedPartOrdsPerProc(numProcs);

  compute_denied_parts_to_remove(ngpMesh, commSparse, localEntities, localRemovePartOrdinals,
                                 deniedKeysPerProc, deniedPartOrdsPerProc);

  stk::pack_and_communicate(partRemovalDenialCommSparse, [&]() {
    for (int destProc = 0; destProc < numProcs; ++destProc) {
      if (destProc == myRank) { continue; }
      auto& sendBuf = partRemovalDenialCommSparse.send_buffer(destProc);

      unsigned numKeysDeniedForRemoval = deniedKeysPerProc[destProc].size();
      sendBuf.pack(numKeysDeniedForRemoval);

      for (unsigned i = 0; i < numKeysDeniedForRemoval; ++i) {
        sendBuf.pack(deniedKeysPerProc[destProc][i]);
        sendBuf.pack(deniedPartOrdsPerProc[destProc][i]);
      }
    }
  });
}

template <typename MeshType>
Kokkos::View<DeniedPartRemoval*, typename MeshType::ngp_mem_space>
unpack_denied_part_removals(MeshType& ngpMesh, CommSparse& partDenialCommSparse)
{
  using EntityKeyWrapper = EntityWrapper<EntityKey>;
  using DeniedViewType = Kokkos::View<DeniedPartRemoval*, typename MeshType::ngp_mem_space>;
  using EntityKeysHostView = Kokkos::View<EntityKeyWrapper*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using PartOrdinalsHostView = Kokkos::View<unsigned*, stk::ngp::HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  auto comm = partDenialCommSparse.parallel();
  auto numProcs = partDenialCommSparse.parallel_size();
  auto myRank = stk::parallel_machine_rank(comm);

  std::vector<EntityKeyWrapper> allDeniedKeys;
  std::vector<unsigned> allDeniedParts;

  for (int srcProc = 0; srcProc < numProcs; ++srcProc) {
    if (srcProc == myRank) { continue; }
    stk::CommBuffer& recvBuf = partDenialCommSparse.recv_buffer(srcProc);
    unsigned numDenied = 0;
    recvBuf.unpack(numDenied);

    for (unsigned i = 0; i < numDenied; ++i) {
      EntityKeyWrapper entityKey;
      unsigned part;
      recvBuf.unpack(entityKey);
      recvBuf.unpack(part);
      allDeniedKeys.push_back(entityKey);
      allDeniedParts.push_back(part);
    }
  }

  if (allDeniedKeys.empty()) {
    return DeniedViewType("deniedPartRemoval", 0);
  }

  EntityKeysHostView deniedKeysHost(allDeniedKeys.data(), allDeniedKeys.size());
  auto deniedKeysDevice = Kokkos::create_mirror_view_and_copy(stk::ngp::MemSpace{}, deniedKeysHost);
  auto deniedEntitiesDevice = ngpMesh.get_entities(deniedKeysDevice);

  PartOrdinalsHostView deniedPartsHost(allDeniedParts.data(), allDeniedParts.size());
  Kokkos::View<unsigned*> deniedPartsDevice("deniedPartsDevice", allDeniedParts.size());
  Kokkos::deep_copy(deniedPartsDevice, deniedPartsHost);

  auto numDeniedResolved = deniedEntitiesDevice.extent(0);
  DeniedViewType deniedPartRemovalView("deniedPartRemoval", numDeniedResolved);
  Kokkos::parallel_for(numDeniedResolved,
    KOKKOS_LAMBDA(const int i) {
      deniedPartRemovalView(i) = DeniedPartRemoval{deniedEntitiesDevice(i), deniedPartsDevice(i)};
    }
  );

  Kokkos::sort(deniedPartRemovalView);

  return deniedPartRemovalView;
}

} } } // namespace stk::mesh::impl

#endif

