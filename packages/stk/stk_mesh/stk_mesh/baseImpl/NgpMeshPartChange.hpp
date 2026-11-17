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

#ifndef STK_MESH_NGPMESHPARTCHANGE_HPP
#define STK_MESH_NGPMESHPARTCHANGE_HPP

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/NgpMeshTypes.hpp"
#include "stk_mesh/baseImpl/NgpViewAlgorithms.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include <cstddef>

namespace stk::mesh::impl {

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

template <typename DeviceBucketRepoType, typename AddPartOrdinalsViewType, typename RemovePartOrdinalsViewType>
bool has_ranked_part(DeviceBucketRepoType const& deviceBucketRepo, AddPartOrdinalsViewType const& addPartOrdinals,
                     RemovePartOrdinalsViewType const& removePartOrdinals)
{
  int numRankedParts = 0;

  Kokkos::parallel_reduce(addPartOrdinals.extent(0),
    KOKKOS_LAMBDA(const int& i, int& update) {
      if (deviceBucketRepo.is_inducible_part(addPartOrdinals(i))) {
        update++;
      }
    }, Kokkos::Sum<int>(numRankedParts)
  );
  Kokkos::fence();

  if (numRankedParts > 0) { return true; }

  Kokkos::parallel_reduce(removePartOrdinals.extent(0),
    KOKKOS_LAMBDA(const int& i, int& update) {
      if (deviceBucketRepo.is_inducible_part(removePartOrdinals(i))) {
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
                              IsLosingPartFn const& isLosingPart)
{
  auto connectedUpperRankEntities = bucket.get_connected_entities(bucketOrd, partRank);
  bool anyKeeps = false;

  Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, connectedUpperRankEntities.size()),
    [&](const int& connIdx, bool& localKeeps) {
      auto upperEntity = connectedUpperRankEntities[connIdx];
      auto bucketId = ngpMesh.device_mesh_index(upperEntity).bucket_id;
      auto& upperBucket = ngpMesh.get_bucket(partRank, bucketId);

      if (!upperBucket.is_owned()) { return; }
      if (!upperBucket.member(partOrdinal)) { return; }
      if (!isLosingPart(upperEntity)) { localKeeps = true; }
    }, Kokkos::LOr<bool>(anyKeeps)
  );

  return anyKeeps;
}


template <typename SortedViewType, typename T>
KOKKOS_INLINE_FUNCTION
bool sorted_view_contains(SortedViewType const& sortedView, T const& target)
{
  const size_t idx = sorted_view_lower_bound(sortedView, target);
  return (idx < sortedView.extent(0)) && (static_cast<T>(sortedView(idx)) == target);
}

template<typename MeshType, typename WrappedEntityViewType, typename AddPartsViewType, typename RmPartsViewType,
         typename NewPartsViewType, typename PartOrdinalsProxyViewType>
void set_new_part_list_per_entity_with_induced_parts(MeshType const& ngpMesh,
                                                     WrappedEntityViewType wrappedEntities,
                                                     AddPartsViewType addParts,
                                                     RmPartsViewType rmParts,
                                                     unsigned maxPartsPerEntity,
                                                     NewPartsViewType newPartsPerEntity,
                                                     PartOrdinalsProxyViewType partOrdinalsProxyView)
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

      auto isInducingPart = [&](Ordinal partOrdinal) { return deviceBucketRepo.is_inducible_part(partOrdinal); };

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
      
      auto isInWrappedEntities = [&](Entity entity) { return sorted_view_contains(wrappedEntities, entity); };

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
              bool keepsPart = should_keep_induced_part(ngpMesh, team, bucket, fastMeshIdx.bucket_ord,
                                                        partToCheck, partToCheckRank, isLosingPart);

              Kokkos::single(Kokkos::PerThread(team), [&]() {
                if (!keepsPart) {
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

// Computes the new part list of each entity in entities and all of its downward-connected entities
template <typename MeshType, typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
auto compute_new_part_lists_with_inducible_parts(MeshType const& ngpMesh,
                                                 const Kokkos::View<EntityWrapper<>*, EntitiesParams...>& entities,
                                                 const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                                 const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
{
  using PartOrdinalsViewType = typename std::remove_reference<decltype(addPartOrdinals)>::type;
  using EntityWrapperViewType = Kokkos::View<EntityWrapper<>*, typename MeshType::ngp_mem_space>;
  using PartOrdinalsProxyViewType = Kokkos::View<PartOrdinalsProxyIndices*>;
  using NewPartListsType = NewPartListsPerEntity<EntityWrapperViewType, PartOrdinalsViewType, PartOrdinalsProxyViewType>;

  Kokkos::Profiling::pushRegion("construct_part_ordinal_proxy");
  Kokkos::Profiling::pushRegion("gather_wrapped_entities");
  auto maxNumDownwardConnectedEntities = get_max_num_downward_connected_entities(ngpMesh, entities);
  auto entityInterval = maxNumDownwardConnectedEntities + 1;
  auto maxNumEntitiesForInducingParts = entityInterval * entities.extent(0);

  EntityWrapperViewType wrappedEntities(Kokkos::view_alloc("wrappedEntities"), maxNumEntitiesForInducingParts);
  populate_all_downward_connected_entities_and_wrap_entities(ngpMesh, entities, entityInterval, wrappedEntities);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("sort_unique_wrapped_entities");
  remove_invalid_wrapped_entities_sort_unique_merge_and_resize(wrappedEntities, typename MeshType::MeshExecSpace{});
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("compute_new_part_lists");
  // determine resulting parts per entity including inducible parts
  const unsigned maxCurrentNumPartsPerEntity = get_max_num_parts_per_entity(ngpMesh, wrappedEntities);
  const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

  PartOrdinalsViewType newPartOrdinalsPerEntity("newPartOrdinals", wrappedEntities.size() * maxNewNumPartsPerEntity);
  PartOrdinalsViewType sortedAddPartOrdinals = get_sorted_view(addPartOrdinals);

  // create part ordinals proxy, sort and unique it (and realloc)
  // A proxy indices view to part ordinals: (rank, startPtr, length)
  PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc("partOrdinalsProxy", Kokkos::WithoutInitializing), wrappedEntities.size());
  set_new_part_list_per_entity_with_induced_parts(ngpMesh, wrappedEntities, sortedAddPartOrdinals, removePartOrdinals,
                                                        maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();

  return NewPartListsType{wrappedEntities, newPartOrdinalsPerEntity, partOrdinalsProxy, maxNewNumPartsPerEntity};
}

// Computes the new part list of each entity in 'entities' for parts that induce nothing. Moves nothing.
template <typename MeshType, typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
auto compute_new_part_lists_without_inducible_parts(MeshType const& ngpMesh,
                                                    const Kokkos::View<EntityWrapper<>*, EntitiesParams...>& entities,
                                                    const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                                    const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
{
  using PartOrdinalsViewType = typename std::remove_reference<decltype(addPartOrdinals)>::type;
  using EntityWrapperViewType = Kokkos::View<EntityWrapper<>*, typename MeshType::ngp_mem_space>;
  using PartOrdinalsProxyViewType = Kokkos::View<PartOrdinalsProxyIndices*>;
  using NewPartListsType = NewPartListsPerEntity<EntityWrapperViewType, PartOrdinalsViewType, PartOrdinalsProxyViewType>;

  Kokkos::Profiling::pushRegion("construct_part_ordinal_proxy");
  const unsigned maxCurrentNumPartsPerEntity = get_max_num_parts_per_entity(ngpMesh, entities);
  const unsigned maxNewNumPartsPerEntity = maxCurrentNumPartsPerEntity + addPartOrdinals.size();

  PartOrdinalsViewType newPartOrdinalsPerEntity("newPartOrdinals", entities.size()*maxNewNumPartsPerEntity);
  PartOrdinalsViewType sortedAddPartOrdinals = get_sorted_view(addPartOrdinals);

  // A proxy indices view to part ordinals: (rank, startPtr, length)
  PartOrdinalsProxyViewType partOrdinalsProxy(Kokkos::view_alloc(Kokkos::WithoutInitializing, "partOrdinalsProxy"), entities.size());
  set_new_part_list_per_entity(ngpMesh, entities, sortedAddPartOrdinals, removePartOrdinals,
                                     maxNewNumPartsPerEntity, newPartOrdinalsPerEntity, partOrdinalsProxy);
  Kokkos::Profiling::popRegion();

  EntityWrapperViewType wrappedEntities = entities;
  return NewPartListsType{wrappedEntities, newPartOrdinalsPerEntity, partOrdinalsProxy, maxNewNumPartsPerEntity};
}

template <typename MeshType, typename... EntitiesParams, typename... AddPartParams, typename... RemovePartParams>
auto compute_local_new_part_lists(MeshType const& ngpMesh,
                                  const Kokkos::View<stk::mesh::Entity*, EntitiesParams...>& entities,
                                  const Kokkos::View<stk::mesh::PartOrdinal*, AddPartParams...>& addPartOrdinals,
                                  const Kokkos::View<stk::mesh::PartOrdinal*, RemovePartParams...>& removePartOrdinals)
{
  Kokkos::Profiling::pushRegion("compute_local_new_part_lists");
  const bool hasRankedPart = has_ranked_part(ngpMesh.get_device_bucket_repository(), addPartOrdinals, removePartOrdinals);
  auto validEntities = remove_invalid_entries_and_resize<typename MeshType::MeshExecSpace>(wrap_entities(entities, true, false),
                                                                                           EntityWrapper(Entity{}));
  auto localNewPartLists = hasRankedPart
                      ? compute_new_part_lists_with_inducible_parts(ngpMesh, validEntities, addPartOrdinals, removePartOrdinals)
                      : compute_new_part_lists_without_inducible_parts(ngpMesh, validEntities, addPartOrdinals, removePartOrdinals);
  Kokkos::Profiling::popRegion();
  return localNewPartLists;
}

} // namespace stk::mesh::impl

#endif
