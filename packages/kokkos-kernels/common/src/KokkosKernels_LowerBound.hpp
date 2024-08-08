//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef _KOKKOSKERNELS_LOWERBOUND_HPP
#define _KOKKOSKERNELS_LOWERBOUND_HPP

/*! \file KokkosKernels_LowerBound.hpp
   Define thread and team-collaborative lower-bound search

   Lower-bound search takes a Kokkos::View, a search value, and a binary
   predicate.
   It returns an index to the first element of the view that does not
   satisfy pred(element, value), or the size of the view if no such
   element exists.

   All elements for which pred(element, value) is true must precede those
   for which it is false.

   The default predicate is less-than, i.e. pred(a,b) = a < b.
   In this case, lower-bound search returns the first index where the value is
   >= the view entry.

   The type of the predicate function must be equivalent to the following:
   \verbatim
   bool operator(const T &a, const T&b);
   \endverbatim
   KokkosKernels_Predicates.hpp defines a variety of common predicates,
   available in KokkosKernels namespace.

   Examples:
   \verbatim
   value  = 3
   view   = {0,1,2,3,4}
          = {t,t,t,f,f}
   result =        3

   value  = -1
   view   = {0,1,2,3,4}
          = {f,f,f,f,f}
   result =  0

   value  = 5
   view   = {0,1,2,3,4}
          = {t,t,t,t,t}
   result =            5

   value  = 1
   view   = {0,1,1,1,2}
          = {t,f,f,f,f}
   result =    1
   \endverbatim

   Contrast with upper-bound, which returns first index for which pred(value,
   element) is true
 */

#include <Kokkos_NumericTraits.hpp>

#include "KokkosKernels_Predicates.hpp"
#include "KokkosKernels_SimpleUtils.hpp"

namespace KokkosKernels {
namespace Impl {

/*! \brief Single-thread sequential lower-bound search

    \tparam ViewLike A Kokkos::View, KokkosKernels::Impl::Iota, or
   KokkosSparse::MergeMatrixDiagonal \tparam Pred a binary predicate function
    \param view the view to search
    \param value the value to search for
    \param pred a binary predicate function
    \returns index of first element in view where pred(element, value) is false,
    or view.size if no such element exists

    At most view.size() predicate function calls
*/
template <typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type>>
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type lower_bound_sequential_thread(
    const ViewLike &view, const typename ViewLike::non_const_value_type &value, Pred pred = Pred()) {
  using size_type = typename ViewLike::size_type;
  static_assert(1 == ViewLike::rank, "lower_bound_sequential_thread requires rank-1 views");

  size_type i = 0;
  while (i < view.size() && pred(view(i), value)) {
    ++i;
  }
  return i;
}

/*! \brief Single-thread binary lower-bound search

    \tparam ViewLike A Kokkos::View, KokkosKernels::Impl::Iota, or
   KokkosSparse::MergeMatrixDiagonal \tparam Pred a binary predicate function
    \param view the view to search
    \param value the value to search for
    \param pred a binary predicate function
    \returns index of first element in view where pred(element, value) is false,
    or view.size if no such element exists

    At most log2(view.size()) + 1 predicate function calls
*/
template <typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type>>
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type lower_bound_binary_thread(
    const ViewLike &view, const typename ViewLike::non_const_value_type &value, Pred pred = Pred()) {
  using size_type = typename ViewLike::size_type;
  static_assert(1 == ViewLike::rank, "lower_bound_binary_thread requires rank-1 views");

  size_type lo = 0;
  size_type hi = view.size();
  while (lo < hi) {
    size_type mid  = (lo + hi) / 2;
    const auto &ve = view(mid);
    if (pred(ve, value)) {  // mid satisfies predicate, look in higher half not
                            // including mid
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

}  // namespace Impl

/*! \brief single-thread lower-bound search

    \tparam ViewLike A Kokkos::View, KokkosKernels::Impl::Iota, or
   KokkosSparse::MergeMatrixDiagonal \tparam Pred a binary predicate function
    \param view the view to search
    \param value the value to search for
    \param pred a binary predicate function
    \returns index of first element in view where pred(element, value) is false,
    or view.size if no such element exists

    This minimizes the calls to predicate:
    for view.size() >= 8, this does a binary search, otherwise, a linear search
*/
template <typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type>>
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type lower_bound_thread(
    const ViewLike &view, const typename ViewLike::non_const_value_type &value, Pred pred = Pred()) {
  static_assert(1 == ViewLike::rank, "lower_bound_thread requires rank-1 views");
  /*
     sequential search makes on average 0.5 * view.size memory accesses
     binary search makes log2(view.size)+1 accesses

     log2(x) <= 0.5x roughly when x >= 8
  */
  if (view.size() >= 8) {
    return Impl::lower_bound_binary_thread(view, value, pred);
  } else {
    return Impl::lower_bound_sequential_thread(view, value, pred);
  }
}

namespace Impl {

/*! \brief Team-collaborative sequential lower-bound search

    \tparam TeamMember the team policy member type
    \tparam ViewLike A Kokkos::View or KokkosKernels::Iota
    \tparam Pred The type of the predicate function to call

    \param handle The Kokkos team handle
    \param view The view-like to search
    \param value The value to compare in the predicate
    \param lo The first index to search
    \param hi One-past the last index to search
    \param pred Apply pred(view(i), value)

    \returns To all team members, the smallest i for which pred(view(i), value)
   is false for i in [lo, hi), or hi if no such value

    Uses a single thread to call \c lower_bound_thread, and broadcasts that
    to all team members.
*/
template <typename TeamMember, typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type>>
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type lower_bound_single_team(
    const TeamMember &handle, const ViewLike &view, const typename ViewLike::non_const_value_type &value,
    Pred pred = Pred()) {
  typename ViewLike::size_type idx;
  Kokkos::single(
      Kokkos::PerTeam(handle),
      [&](typename ViewLike::size_type &lidx) { lidx = KokkosKernels::lower_bound_thread(view, value, pred); }, idx);
  return idx;
}

/*! \brief Team-collaborative sequential lower-bound search

    \tparam TeamMember the team policy member type
    \tparam ViewLike A Kokkos::View or KokkosKernels::Iota
    \tparam Pred The type of the predicate function to call

    \param handle The Kokkos team handle
    \param view The view-like to search
    \param value The value to compare in the predicate
    \param lo The first index to search
    \param hi One-past the last index to search
    \param pred Apply pred(view(i), value)

    \returns To all team members, the smallest i for which pred(view(i), value)
   is false for i in [lo, hi), or hi if no such value

    Apply pred(view(i), value) for i in [lo, hi)
*/
template <typename TeamMember, typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type>>
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type lower_bound_sequential_team(
    const TeamMember &handle, const ViewLike &view, const typename ViewLike::non_const_value_type &value,
    typename ViewLike::size_type lo, typename ViewLike::size_type hi, Pred pred = Pred()) {
  using size_type = typename ViewLike::size_type;
  static_assert(1 == ViewLike::rank, "lower_bound_sequential_team requires rank-1 views");
  static_assert(is_iota_v<ViewLike> || Kokkos::is_view<ViewLike>::value,
                "lower_bound_sequential_team requires a "
                "KokkosKernels::Impl::Iota or a Kokkos::View");

  if (lo == hi) {
    return hi;
  }
  size_type teamI;
  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(handle, lo, hi),
      [&](const size_type &i, size_type &li) {
        li = KOKKOSKERNELS_MACRO_MIN(li, hi);
        if (i < li) {                   // no need to search higher than the smallest so far
          if (!pred(view(i), value)) {  // look for the smallest index that does
                                        // not satisfy
            li = i;
          }
        }
      },
      Kokkos::Min<size_type>(teamI));
  return teamI;
}

/*! \brief Team-collaborative sequential lower-bound search

    \tparam TeamMember the team policy member type
    \tparam ViewLike A Kokkos::View or KokkosKernels::Iota
    \tparam Pred The type of the predicate function to call

    \param handle The Kokkos team handle
    \param view The view-like to search
    \param value The value to compare in the predicate
    \param pred Apply pred(view(i), value)

    \returns To all team members, the smallest i for which pred(view(i), value)
   is false or view.size() if no such value
*/
template <typename TeamMember, typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type>>
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type lower_bound_sequential_team(
    const TeamMember &handle, const ViewLike &view, const typename ViewLike::non_const_value_type &value,
    Pred pred = Pred()) {
  return lower_bound_sequential_team(handle, view, value, 0, view.size(), pred);
}

/*! \brief A range for the k-ary lower bound search

    The RangeReducer will maximize the lower bound and
    minimize the upper bound
*/
template <typename T>
struct Range {
  T lb;  /// lower-bound
  T ub;  /// upper-bound

  KOKKOS_INLINE_FUNCTION
  Range() { init(); }

  KOKKOS_INLINE_FUNCTION
  constexpr Range(const T &_lb, const T &_ub) : lb(_lb), ub(_ub) {}

  KOKKOS_INLINE_FUNCTION
  void init() {
    lb = Kokkos::Experimental::finite_min_v<T>;  // will be max'd
    ub = Kokkos::Experimental::finite_max_v<T>;  // will be min'd
  }
};

/// \brief maximizes the lower bound, and minimizes the upper bound of a Range
template <typename T, typename Space>
struct RangeReducer {
  using reducer          = RangeReducer;
  using value_type       = Range<T>;
  using result_view_type = Kokkos::View<Range<T> *, Space, Kokkos::MemoryUnmanaged>;

 private:
  value_type &value;

 public:
  KOKKOS_INLINE_FUNCTION
  RangeReducer(value_type &value_) : value(value_) {}

  KOKKOS_INLINE_FUNCTION
  void join(value_type &dst, const value_type &src) const {
    dst.lb = KOKKOSKERNELS_MACRO_MAX(dst.lb, src.lb);
    dst.ub = KOKKOSKERNELS_MACRO_MIN(dst.ub, src.ub);
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &val) const { val.init(); }

  KOKKOS_INLINE_FUNCTION
  value_type &reference() const { return value; }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return result_view_type(&value, 1); }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return true; }
};

/*! \brief team-collaborative K-ary lower-bound search

    \tparam TeamMember the team policy member type
    \tparam ViewLike A Kokkos::View or KokkosKernels::Iota
    \tparam Pred the binary predicate function type

    Actually, K+1-ary, where K is the size of the team
    Split the view into k+1 segments at K points
    Evalute the predicate in parallel at each point and use a joint min-max
   parallel reduction:
      * The lower bound is after the max index where the predicate was true
      * The upper bound is no greater than the min index where the predicate was
   false Once there are fewer values left than threads in the team, switch to
   team sequential search
*/
template <typename TeamMember, typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type>>
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type lower_bound_kary_team(
    const TeamMember &handle, const ViewLike &view, const typename ViewLike::non_const_value_type &value,
    Pred pred = Pred()) {
  static_assert(1 == ViewLike::rank, "lower_bound_kary_team requires rank-1 views");
  static_assert(is_iota_v<ViewLike> || Kokkos::is_view<ViewLike>::value,
                "lower_bound_kary_team requires a "
                "KokkosKernels::Impl::Iota or a Kokkos::View");

  using size_type = typename ViewLike::size_type;

  size_type lo = 0;
  size_type hi = view.size();
  while (lo < hi) {
    // if fewer than team_size elements left, just hit them all sequentially
    if (lo + handle.team_size() >= hi) {
      return lower_bound_sequential_team(handle, view, value, lo, hi, pred);
    }

    // otherwise, split the region up among threads
    size_type mid = lo + (hi - lo) * (handle.team_rank() + 1) / (handle.team_size() + 1);
    auto ve       = view(mid);

    // reduce across threads to figure out where the new search bounds are
    // if a thread satisfies the predicate, the first element that does not
    // satisfy must be after that thread's search point. we want the max such
    // point across all threads if a thread does not satisfy the predicate, the
    // first element that does not satisfy must be before or equal. we want the
    // min such point across all threads
    Range<size_type> teamRange;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(handle, 0, handle.team_size()),
        [&](const int &, Range<size_type> &lr) {
          lr.lb = KOKKOSKERNELS_MACRO_MAX(lo, lr.lb);  // no lower than lo
          lr.ub = KOKKOSKERNELS_MACRO_MIN(hi, lr.ub);  // no higher than hi
          // if pred(view(mid), value), then the lower bound is above this
          if (pred(ve, value)) {
            lr.lb = mid + 1;
          } else {  // otherwise the lower bound is no larger than this
            lr.ub = mid;
          }
        },
        RangeReducer<size_type, typename ViewLike::device_type>(teamRange));

    // next iteration, search in the newly-discovered window
    hi = teamRange.ub;
    lo = teamRange.lb;
  }
  return lo;
}

}  // namespace Impl

/*! \brief Team-collaborative lower-bound search

    \tparam TeamMember the team policy member type the Kokkos team handle
    \tparam View the type of view
    \tparam Pred the type of the predicate

    \param handle a Kokkos team handle
    \param view a Kokkos::View to search
    \param value the value to search for
    \param pred the predicate to test entries in the view

    \returns The smallest i in range [0, view.size()) for which pred(view(i),
   value) is not true, or view.size() if no such `i` exists

    default pred is `element < value`, i.e. return the index to the first
   element in the view that does not satisfy `element < value`. For well-ordered
   types this is the first element where element >= value

    Pred should be a binary function comparing two `typename
   View::non_const_value_type`
*/
template <typename TeamMember, typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type>>
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type lower_bound_team(
    const TeamMember &handle, const ViewLike &view, const typename ViewLike::non_const_value_type &value,
    Pred pred = Pred()) {
  static_assert(1 == ViewLike::rank, "lower_bound_team requires rank-1 views");
  static_assert(KokkosKernels::Impl::is_iota_v<ViewLike> || Kokkos::is_view<ViewLike>::value,
                "lower_bound_team requires a "
                "KokkosKernels::Impl::Iota or a Kokkos::View");

  /* kary search is A = (k-1) * (logk(view.size()) + 1) accesses

     sequential search is B = view.size() accesses

      A < B is true ruoughly when view.size() > 3 * k
  */
  if (view.size() > 3 * size_t(handle.team_size())) {
    return Impl::lower_bound_kary_team(handle, view, value, pred);
  } else {
    return Impl::lower_bound_sequential_team(handle, view, value, pred);
  }
}

}  // namespace KokkosKernels

#endif  // _KOKKOSKERNELS_LOWERBOUND_HPP