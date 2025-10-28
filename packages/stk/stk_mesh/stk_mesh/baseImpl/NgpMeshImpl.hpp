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

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/baseImpl/ViewVector.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include "Kokkos_Core.hpp"
#include "Kokkos_Macros.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

namespace stk {
namespace mesh {
namespace impl {

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
      if (lhs(i) < rhs(i)) {   // FIXME
        return true;
      }
    }
    return false;
  }
};

struct EntityWrapper
{
  KOKKOS_FUNCTION
  EntityWrapper()
    : entity(Entity{}),
      isUserInputEntity(false),
      isForPartInduction(false)
  {}

  KOKKOS_FUNCTION
  EntityWrapper(Entity e, bool u = false, bool p = false)
    : entity(e),
      isUserInputEntity(u),
      isForPartInduction(p)
  {}

  KOKKOS_INLINE_FUNCTION
  operator Entity() const {
    return entity;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(EntityWrapper other) const {
    return (entity == other.entity);
  }

  KOKKOS_INLINE_FUNCTION
  bool operator<(EntityWrapper other) const {
    if (entity == other.entity) {
      return isForPartInduction == false;
    } else {
      return entity < other;
    }
  }

  Entity entity;
  bool isUserInputEntity;
  bool isForPartInduction;
};

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

template <typename T>
T get_round_up_pow_of_2(T val)
{
  static_assert(std::is_integral_v<T>);

  if (val <= 0) { return val; }

  val--;
  val |= val >> 1;
  val |= val >> 2;
  val |= val >> 4;
  val |= val >> 8;
  val |= val >> 16;
  val++;

  return val;
}

template <typename ViewType, typename ExecSpace>
void sort_and_unique(ViewType& view, ExecSpace const& execSpace)
{
  Kokkos::sort(view);
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, view));
  Kokkos::Experimental::unique(execSpace, view);
}

template <typename ViewType, typename ExecSpace>
void sort_and_unique_and_resize(ViewType& view, ExecSpace const& execSpace)
{
  Kokkos::sort(view);
  STK_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, view));

  auto newEnd = Kokkos::Experimental::unique(execSpace, view);
  auto begin = Kokkos::Experimental::begin(view);
  size_t newSize = Kokkos::Experimental::distance(begin, newEnd);
  if (newSize != view.extent(0)) {
    Kokkos::resize(view, newSize);
  }
}

template <typename DeviceBucketRepoType, typename PartOrdinalsViewType>
bool has_ranked_part(DeviceBucketRepoType const& deviceBucketRepo, PartOrdinalsViewType const& addPartOrdinals, PartOrdinalsViewType const& removePartOrdinals)
{
  int numRankedParts = 0;

  Kokkos::parallel_reduce(addPartOrdinals.extent(0),
    KOKKOS_LAMBDA(const int& i, int& update) {
      if (deviceBucketRepo.is_ranked_part(addPartOrdinals(i))) {
        update++;
      }
    }, Kokkos::Sum<int>(numRankedParts)
  );
  Kokkos::fence();

  if (numRankedParts > 0) { return true; }

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
        totalNumConnectedEntities += bucket.get_connected_entities(0, rankToSearch).size();
      }

      update = (update > totalNumConnectedEntities) ? update : totalNumConnectedEntities;
    }, Kokkos::Max<int>(max)
  );
  Kokkos::fence();

  return max;
}

template <typename MeshType, typename EntityViewType, typename WrappedEntityViewType>
void populate_all_downward_connected_entities_and_wrap_entities(MeshType const& ngpMesh, EntityViewType const& entities,
                                                                int entityInterval, WrappedEntityViewType const& wrappedEntities)
{
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int i) {
      auto myStartIdx = i * entityInterval;
      auto myCurrentIdx = myStartIdx;
      auto rank = ngpMesh.entity_rank(entities(i));
      auto fastMeshIdx = ngpMesh.device_mesh_index(entities(i));
      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIdx.bucket_id);

      wrappedEntities(myCurrentIdx++) = EntityWrapper(entities(i), true);

      for (auto iRank = rank; iRank >= stk::topology::NODE_RANK; --iRank) {
        if (iRank == rank) { continue; }
        auto connectedEntities = bucket.get_connected_entities(fastMeshIdx.bucket_ord, iRank);

        for (unsigned j = 0; j < connectedEntities.size(); j++) {
          auto wrappedEntity = EntityWrapper(connectedEntities[j], false, true);
          wrappedEntities(myCurrentIdx++) = wrappedEntity;
        }
      }
    }
  );
  Kokkos::fence();
}

template <typename EntityViewType, typename ExecSpace>
void remove_invalid_entities_sort_unique_and_resize(EntityViewType& wrappedEntities, ExecSpace const& execSpace)
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

  auto newEnd = Kokkos::Experimental::unique(execSpace, uview);
  auto newBegin = Kokkos::Experimental::begin(uview);
  size_t newSize = Kokkos::Experimental::distance(newBegin, newEnd);
  if (newSize != wrappedEntities.extent(0)) {
    Kokkos::resize(wrappedEntities, newSize);
  }
}

template<typename MeshType, typename WrappedEntityViewType, typename AddPartsViewType, typename RmPartsViewType,
         typename NewPartsViewType, typename PartOrdinalsProxyViewType>
void set_new_part_list_per_entity_with_induced_parts(const MeshType& ngpMesh,
                                                     const WrappedEntityViewType& wrappedEntities,
                                                     const AddPartsViewType& addParts,
                                                     const RmPartsViewType& rmParts,
                                                     unsigned maxPartsPerEntity,
                                                     NewPartsViewType& newPartsPerEntity,
                                                     PartOrdinalsProxyViewType& partOrdinalsProxyView)
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

      auto isRemovedPart = [&](Ordinal partOrdinal) { return partOrdinal == InvalidPartOrdinal; };

      auto isRankedPart = [&](Ordinal partOrdinal) { return deviceBucketRepo.is_ranked_part(partOrdinal); };

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
          auto rankedPart = isRankedPart(partToCheck);

          if (!rankedPart && inAddParts) {
            if (!isUserInputEntity) {
              newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
              return;
            }
          } else if(!rankedPart && inRmParts) {
            if (isUserInputEntity) {
              newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
              return;
            }
          }
          else if (rankedPart && inAddParts) {
            if (rank == partToCheckRank) {
              return;
            } else if (partToCheckRank > rank) {
              if (!isForPartInduction) {
                newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
                return;
              }
            }
          } else if (rankedPart && inRmParts) {
            if (rank == partToCheckRank) {
              newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
              return;
            } else if (partToCheckRank > rank) {
              auto connectedUpperRankEntities = bucket.get_connected_entities(0, partToCheckRank);
              bool foundAllInWrappedEntities = true;

              Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, connectedUpperRankEntities.size()),
                [&](const int& connIdx, bool& localFound) {
                  auto entity = connectedUpperRankEntities[connIdx];
                  auto bucketId = ngpMesh.device_mesh_index(entity).bucket_id;
                  auto& localBucket = ngpMesh.get_bucket(partToCheckRank, bucketId);
                  auto isMember = localBucket.member(partToCheck);

                  if (!isMember) {
                    localFound &= true;
                    return;
                  }

                  auto found = isInWrappedEntities(entity);
                  localFound &= found;
                }, Kokkos::LAnd<bool>(foundAllInWrappedEntities)
              );

              Kokkos::single(Kokkos::PerThread(team), [&]() {
                if (foundAllInWrappedEntities) {
                  newPartsPerEntity(myPartOrdinalIdx) = InvalidPartOrdinal;
                }
              });
              return;
            }
          }
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
} } } // namespace stk::mesh::impl

#endif

