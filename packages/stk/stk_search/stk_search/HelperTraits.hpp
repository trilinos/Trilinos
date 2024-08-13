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

#ifndef STK_SEARCH_HELPER_TRAITS_H
#define STK_SEARCH_HELPER_TRAITS_H

#include <type_traits>
#include "BoxIdent.hpp"
#include "Kokkos_Core.hpp"
#include "stk_search/SearchMethod.hpp"

#include <vector>

namespace stk::search {

template <typename T>
struct is_box_ident : std::false_type
{};

template <typename BoxType, typename IdentType>
struct is_box_ident<BoxIdent<BoxType, IdentType>> : std::true_type
{};

template <typename T>
struct is_box_ident_proc : std::false_type
{};

template <typename BoxType, typename IdentProcType>
struct is_box_ident_proc<BoxIdentProc<BoxType, IdentProcType>> : std::true_type
{};

template <typename T>
struct is_ident_intersection : std::false_type
{};

template <typename DomainIdent, typename RangeIdent>
struct is_ident_intersection<IdentIntersection<DomainIdent, RangeIdent>> : std::true_type
{};

template <typename T>
struct is_ident_proc_intersection : std::false_type
{};

template <typename DomainIdentProc, typename RangeIdentProc>
struct is_ident_proc_intersection<IdentProcIntersection<DomainIdentProc, RangeIdentProc>> : std::true_type
{};


template <typename T, typename = void>
struct value_type_or_void
{
  using type = void;
};

template <typename T>
struct value_type_or_void<T, std::void_t<typename T::value_type>>
{
  using type = typename T::value_type;
};

template <typename T>
using value_type_or_void_t = typename value_type_or_void<std::remove_reference_t<T>>::type;


// replace with std::remove_cvref in c++20
template <typename T>
struct remove_cvref
{
  using type = std::remove_cv_t<std::remove_reference_t<T>>;
};

template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

template <typename T>
constexpr bool is_box_ident_v = is_box_ident<remove_cvref_t<T>>::value;

template <typename T>
constexpr bool is_box_ident_container_v = is_box_ident_v<value_type_or_void_t<remove_cvref_t<T>>>;

template <typename T>
constexpr bool is_box_ident_proc_v = is_box_ident_proc<remove_cvref_t<T>>::value;

template <typename T>
constexpr bool is_box_ident_proc_container_v = is_box_ident_proc_v<value_type_or_void_t<remove_cvref_t<T>>>;

template <typename T>
constexpr bool is_ident_intersection_v = is_ident_intersection<remove_cvref_t<T>>::value;

template <typename T>
constexpr bool is_ident_intersection_container_v = is_ident_intersection_v<value_type_or_void_t<T>>;

template <typename T>
constexpr bool is_ident_proc_intersection_v = is_ident_proc_intersection<remove_cvref_t<T>>::value;

template <typename T>
constexpr bool is_ident_proc_intersection_container_v = is_ident_proc_intersection_v<value_type_or_void_t<T>>;

template <typename T>
constexpr bool is_modifiable_v = !std::is_const_v<std::remove_reference_t<T>>;

template <typename T>
constexpr bool is_modifiable_view_v = Kokkos::is_view_v<typename std::remove_reference_t<T>> &&
                                      is_modifiable_v<typename std::remove_reference_t<T>::value_type>;

template <typename ExecSpaceT, typename ViewTypeT>
constexpr void check_view_is_usable_from()
{
  using ExecSpace = std::remove_cv_t<ExecSpaceT>;
  using ViewType = std::remove_cv_t<ViewTypeT>;
  static_assert(Kokkos::is_execution_space_v<ExecSpace>, "type passed as ExecSpace is not an ExecutionSpace");
  static_assert(Kokkos::is_view_v<ViewType>, "type passed in as ViewType is not a view");
  static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename ViewType::memory_space>::accessible, "ViewType is not accesible from ExecSpace");
}

template <typename DomainOrRangeView, typename ExecutionSpace>
constexpr void check_domain_or_range_view_types_local()
{
  check_view_is_usable_from<ExecutionSpace, DomainOrRangeView>();
  static_assert(DomainOrRangeView::rank() == 1, "View must be a rank 1 View");
  static_assert(is_box_ident_container_v<DomainOrRangeView>, "View must be a View of BoxIdent<T, U>");
}


template <typename DomainView, typename RangeView, typename ResultView, typename ExecutionSpace>
constexpr void check_coarse_search_types_local()
{
  check_domain_or_range_view_types_local<DomainView, ExecutionSpace>();
  check_domain_or_range_view_types_local<RangeView, ExecutionSpace>();
  check_view_is_usable_from<ExecutionSpace, ResultView>();
  static_assert(ResultView::rank() == 1, "ResultView must be a rank 1 View");
  static_assert(is_ident_intersection_container_v<ResultView>, "ResultView must be a View of IdentIntersection<T, U>");
  static_assert(is_modifiable_view_v<ResultView>, "ResultView must not be const (ie. View<IdentProcIntersections> not View<const IdentProcIntersections>)");
}

template <typename DomainOrRangeView, typename ExecutionSpace>
constexpr void check_domain_or_range_view_parallel()
{
  check_view_is_usable_from<ExecutionSpace, DomainOrRangeView>();
  static_assert(DomainOrRangeView::rank() == 1, "View must be a rank 1 View");
  static_assert(is_box_ident_proc_container_v<DomainOrRangeView>, "View must be a View of BoxIdentProc<T, U>");
}

template <typename DomainView, typename RangeView, typename ResultView, typename ExecutionSpace>
constexpr void check_coarse_search_types_parallel()
{
  check_domain_or_range_view_parallel<DomainView, ExecutionSpace>();
  check_domain_or_range_view_parallel<RangeView, ExecutionSpace>();
  check_view_is_usable_from<ExecutionSpace, ResultView>();
  static_assert(ResultView::rank() == 1, "ResultView must be a rank 1 View");
  static_assert(is_ident_proc_intersection_container_v<ResultView>, "ResultView must be a View of IdentProcIntersection<T, U>");
  static_assert(is_modifiable_view_v<ResultView>, "ResultView must not be const (ie. View<IdentProcIntersections> not View<const IdentProcIntersections>)");
}

}

#endif