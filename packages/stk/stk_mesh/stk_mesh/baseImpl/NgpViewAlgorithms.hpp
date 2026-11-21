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

#ifndef STK_MESH_NGPVIEWALGORITHMS_HPP
#define STK_MESH_NGPVIEWALGORITHMS_HPP

#include "stk_mesh/base/Entity.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include <cstddef>

namespace stk::mesh::impl {

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

template <typename SortedViewType, typename T>
KOKKOS_INLINE_FUNCTION
size_t sorted_view_lower_bound(SortedViewType const& sortedView, T const& target)
{
  size_t lo = 0;
  size_t hi = sortedView.extent(0);
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    if (static_cast<T>(sortedView(mid)) < target) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

} // namespace stk::mesh::impl

#endif
