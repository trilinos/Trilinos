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
#include <Kokkos_Sort.hpp>
#include <Kokkos_StdAlgorithms.hpp>

namespace stk {
namespace mesh {
namespace impl {

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
      lmax = partOrdinalsPair.second - partOrdinalsPair.first;
    }, Kokkos::Max<unsigned>(max)
  );

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

template<typename MESH_TYPE, typename EntityViewType, typename AddPartsViewType, typename RmPartsViewType, typename NewPartsViewType>
void set_new_part_list_per_entity(const MESH_TYPE& ngpMesh, const EntityViewType& entities,
                                  const AddPartsViewType& addParts,
                                  const RmPartsViewType& rmParts,
                                  unsigned maxPartsPerEntity,
                                  NewPartsViewType& newPartsPerEntity)
{
  using BucketType = typename MESH_TYPE::BucketType;
  using ExecSpace = typename MESH_TYPE::MeshExecSpace;
  using TeamMember = typename stk::ngp::TeamPolicy<ExecSpace>::member_type;

  STK_NGP_ThrowAssert(Kokkos::Experimental::is_sorted(ExecSpace{}, addParts));
  auto teamPolicy = stk::ngp::TeamPolicy<ExecSpace>(entities.size(), Kokkos::AUTO);

  Kokkos::parallel_for("set_new_part_lists", teamPolicy,
    KOKKOS_LAMBDA(const TeamMember& teamMember) {
      const unsigned i = teamMember.league_rank();
      const unsigned myStartIdx = i*(maxPartsPerEntity+1);
      newPartsPerEntity(myStartIdx) = 0;
      const stk::mesh::EntityRank rank = ngpMesh.entity_rank(entities(i));
      const stk::mesh::FastMeshIndex entityIdx = ngpMesh.device_mesh_index(entities(i));
      const BucketType& bucket = ngpMesh.get_bucket(rank, entityIdx.bucket_id);
      auto currentPartOrdsPair = bucket.superset_part_ordinals();

      const unsigned idx = myStartIdx+1;
      Ordinal* dest = &newPartsPerEntity(idx);
      const Ordinal* first1 = addParts.data();
      const Ordinal* last1 = first1+addParts.size();
      const Ordinal* first2 = currentPartOrdsPair.first;
      const Ordinal* last2 = currentPartOrdsPair.second;
      unsigned myNumParts = (last2-first2) + addParts.size();
      Kokkos::single(Kokkos::PerTeam(teamMember),[&]() {
        my_merge(first1, last1, first2, last2, dest);
        newPartsPerEntity(myStartIdx) = myNumParts;
      });

      auto isInRmParts = [&](Ordinal item) {
                           for(unsigned rp=0; rp<rmParts.size(); ++rp)  {
                             if (item == rmParts(rp)) { return true; }
                           }
                           return false;
                         };

      auto begin = Kokkos::Experimental::begin(newPartsPerEntity)+idx;
      auto end = begin + myNumParts;

      end = Kokkos::Experimental::remove_if(teamMember, begin, end,isInRmParts);
      end = Kokkos::Experimental::unique(teamMember, begin, end);

      Kokkos::single(Kokkos::PerTeam(teamMember),[&]() {
        myNumParts = end - begin;
        newPartsPerEntity(myStartIdx) = myNumParts;
      });
    }
  );
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

} } } // namespace stk::mesh::impl

#endif

