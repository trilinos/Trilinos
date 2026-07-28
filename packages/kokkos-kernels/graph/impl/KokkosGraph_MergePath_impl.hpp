// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef _KOKKOSGRAPH_MERGEPATH_IMPL_HPP
#define _KOKKOSGRAPH_MERGEPATH_IMPL_HPP

#include "KokkosKernels_Iota.hpp"
#include "KokkosKernels_LowerBound.hpp"

/*! \file KokkosGraph_MergePath_impl.hpp
 *
 */

namespace KokkosGraph {
namespace Impl {

/* the index in a and b in diagonal_search where the diagonal crosses the merge
 * path
 */
template <typename a_index_type, typename b_index_type>
struct DiagonalSearchResult {
  KOKKOS_INLINE_FUNCTION DiagonalSearchResult() : ai(0), bi(0) {}
  a_index_type ai;
  b_index_type bi;
};

/*! \brief a view into the entries of the Merge Matrix along a diagonal

   \tparam AView: A Kokkos::View, the type of A
   \tparam BViewLike A Kokkos::View or an Iota, the type of B

   Example merge matrix M of two arrays
   A (vertical) and B (horizontal),
   as seen in Odeh, Green, Mwassi, Shmueli, Birk
   Merge Path - Parallel Merging Made Simple
   2012

   M[i,j] = 1 iff A[i] > B[j]

   operator(k) returns A[i] > B[j] at the kth entry of the diagonal

        3  5 12 22 45 64 69 82
      ------------------------
      |  /           /
   17 | 1  1  1  0  0  0  0  0
      |/          /
   29 | 1  1  1  1  0  0  0  0
      |         /
   35 | 1  1  1  1  0  0  0  0
      |     /
   73 | 1  1  1  1  1  1  1  0
      |   /
   86 | 1  1  1  1  1  1  1  1
      |/
   90 | 1  1  1  1  1  1  1  1
      |
   95 | 1  1  1  1  1  1  1  1
      |
   99 | 1  1  1  1  1  1  1  1

  Diagonals are counted from the top-left.
  Diagonals are indexed from the bottom-left.
  Shown on the figure above is the 1st and 5th diagonal

  The 0th diagonal D_0 has length 0
  The 1st diagonal D_1 has length 1
  The 5th diagonal D_5 has length 5
  The 9th diagonal D_9 has length 7

  D_1(0) = 1
  D_5(0..3) = 1
  D_5(4) = 0
*/
template <typename AView, typename BViewLike>
class MergeMatrixDiagonal {
 public:
  static_assert(AView::rank == 1, "MergeMatrixDiagonal AView must be rank 1");
  static_assert(BViewLike::rank == 1, "MergeMatrixDiagonal BViewLike must be rank 1");

  // implement bare minimum parts of the view interface
  enum { rank = 1 };
  typedef bool non_const_value_type;
  typedef typename std::conditional<
      sizeof(typename AView::size_type) >= sizeof(typename BViewLike::size_type), typename AView::size_type,
      typename BViewLike::size_type>::type size_type;  // larger size_type of the two view-like types

  using result_type = DiagonalSearchResult<typename AView::size_type, typename BViewLike::size_type>;

  KOKKOS_INLINE_FUNCTION
  MergeMatrixDiagonal(const AView &a, const BViewLike &b, const size_type diagonal) : a_(a), b_(b), d_(diagonal) {}
  MergeMatrixDiagonal() = default;

  KOKKOS_INLINE_FUNCTION
  result_type result(const size_type &di) const noexcept {
    result_type res;
    if (0 == d_) {
      res.ai = 0;
      res.bi = 0;
      return res;
    } else {
      res = diag_to_a_b(di);
      res.ai += 1;
      return res;
    }
  }

  // compare a[i] > b[j] along diagonal at entry di
  KOKKOS_INLINE_FUNCTION
  bool operator()(const size_type di) const {
    result_type res = diag_to_a_b(di);
    if (size_t(res.ai) >= a_.size()) {
      return true;  // on the +a side out of matrix bounds is 1
    } else if (size_t(res.bi) >= b_.size()) {
      return false;  // on the +b side out of matrix bounds is 0
    } else {
      return a_(res.ai) > b_(res.bi);
    }
  }

  /*! \brief length of the diagonal

  */
  KOKKOS_INLINE_FUNCTION
  size_type size() const noexcept {
    if (d_ <= a_.size() && d_ <= b_.size()) {
      return d_;
    } else if (d_ > a_.size() && d_ > b_.size()) {
      // TODO: this returns nonsense if d_ happens to be outside the merge
      // matrix
      return a_.size() + b_.size() - d_;
    } else {
      return KOKKOSKERNELS_MACRO_MIN(a_.size(), b_.size());
    }
  }

 private:
  // translate an index along the diagonal to indices into a_ and b_
  KOKKOS_INLINE_FUNCTION
  result_type diag_to_a_b(const size_type &di) const noexcept {
    result_type res;
    res.ai = d_ < a_.size() ? (d_ - 1) - di : a_.size() - 1 - di;
    res.bi = d_ < a_.size() ? di : d_ + di - a_.size();
    return res;
  }

  AView a_;
  BViewLike b_;
  size_type d_;  // diagonal
};

/*! \brief Return the first index on diagonal \code diag
           in the merge matrix of \code a and \code b that is not 1

This is effectively a lower-bound search on the merge matrix diagonal
where the predicate is "equals 1"
*/
template <typename AView, typename BViewLike>
KOKKOS_INLINE_FUNCTION DiagonalSearchResult<typename AView::size_type, typename BViewLike::size_type> diagonal_search(
    const AView &a, const BViewLike &b, typename MergeMatrixDiagonal<AView, BViewLike>::size_type diag) {
  // unmanaged view types for a and b
  typedef Kokkos::View<typename AView::value_type *, typename AView::device_type, Kokkos::MemoryUnmanaged> um_a_view;
  typedef Kokkos::View<typename BViewLike::value_type *, typename BViewLike::device_type, Kokkos::MemoryUnmanaged>
      um_b_view;

  um_a_view ua(a.data(), a.size());

  // if BViewLike is an Iota, pass it on directly to MMD,
  // otherwise, create an unmanaged view of B
  typedef typename std::conditional<KokkosKernels::Impl::is_iota<BViewLike>::value, BViewLike, um_b_view>::type b_type;

  typedef MergeMatrixDiagonal<um_a_view, b_type> MMD;
  MMD mmd;
  if constexpr (KokkosKernels::Impl::is_iota<BViewLike>::value) {
    mmd = MMD(ua, b, diag);
  } else {
    b_type ub(b.data(), b.size());
    mmd = MMD(ua, ub, diag);
  }

  // returns index of the first element that does not satisfy pred(element,
  // value) our input view is the merge matrix entry along the diagonal, and we
  // want the first one that is not true. so our predicate just tells us if the
  // merge matrix diagonal entry is equal to true or not
  const typename MMD::size_type idx = KokkosKernels::lower_bound_thread(mmd, true, KokkosKernels::Equal<bool>());
  return mmd.result(idx);
}

template <typename TeamMember, typename AView, typename BViewLike>
KOKKOS_INLINE_FUNCTION DiagonalSearchResult<typename AView::size_type, typename BViewLike::size_type> diagonal_search(
    const TeamMember &handle, const AView &a, const BViewLike &b,
    typename MergeMatrixDiagonal<AView, BViewLike>::size_type diag) {
  // unmanaged view types for a and b
  typedef Kokkos::View<typename AView::value_type *, typename AView::device_type, Kokkos::MemoryUnmanaged> um_a_view;
  typedef Kokkos::View<typename BViewLike::value_type *, typename BViewLike::device_type, Kokkos::MemoryUnmanaged>
      um_b_view;

  um_a_view ua(a.data(), a.size());

  // if BViewLike is an Iota, pass it on directly to MMD,
  // otherwise, create an unmanaged view of B
  typedef typename std::conditional<KokkosKernels::Impl::is_iota<BViewLike>::value, BViewLike, um_b_view>::type b_type;

  typedef MergeMatrixDiagonal<um_a_view, b_type> MMD;
  MMD mmd;
  if constexpr (KokkosKernels::Impl::is_iota<BViewLike>::value) {
    mmd = MMD(ua, b, diag);
  } else {
    b_type ub(b.data(), b.size());
    mmd = MMD(ua, ub, diag);
  }

  // returns index of the first element that does not satisfy pred(element,
  // value) our input view is the merge matrix entry along the diagonal, and we
  // want the first one that is not true. so our predicate just tells us if the
  // merge matrix diagonal entry is equal to true or not
  const typename MMD::size_type idx = KokkosKernels::lower_bound_team(handle, mmd, true, KokkosKernels::Equal<bool>());
  return mmd.result(idx);
}

/* search for the intersection of diagonal d into a and b

   a and b are allowed to be different types so one of them can be an Iota

   returns DiagonalSearchResult where `ai` is the index in a and `bi` is the
   index in `b` where the diagonal crosses the merge path

   If the diagonal requested is too large, `ai` = a.size() and `bi` = b.size()
   is returned
*/
template <typename AView, typename BView>
KOKKOS_INLINE_FUNCTION DiagonalSearchResult<typename AView::size_type, typename BView::size_type> diagonal_search2(
    const AView &a, const BView &b, typename AView::size_type diag) {
  typedef typename AView::size_type size_type;
  typedef typename AView::non_const_value_type ordinal_t;
  DiagonalSearchResult<typename AView::size_type, typename BView::size_type> res;

  if (diag >= a.size() + b.size()) {  // diagonal outside the grid
    res.ai = a.size();
    res.bi = b.size();
    return res;
  }

  size_type lo = 0;
  // hi - lo is the length of diagonal `diag`
  size_type hi;
  if (diag <= a.size() && diag <= b.size()) {
    hi = diag;
  } else if (diag > a.size() && diag > b.size()) {
    hi = a.size() + b.size() - diag;
  } else {
    hi = KOKKOSKERNELS_MACRO_MIN(a.size(), b.size());
  }

  // fprintf(stderr, "lo=%ld hi=%ld\n", lo, hi);
  // diagonal is indexed in the positive b and negative a direction
  while (hi > lo) {
    size_type mid = (lo + hi) / 2;
    size_type ai  = diag <= a.size() ? diag - mid - 1 : a.size() - mid - 1;
    size_type bi  = diag <= a.size() ? mid : diag - a.size() + mid;

    // printf("lo=%ld hi=%ld mid=%ld ai=%ld bi=%ld\n", lo, hi, mid, ai, bi);

    const ordinal_t av = a(ai);
    const ordinal_t bv = b(bi);
    // std::cerr << "av=" << av << " bv=" << bv << std::endl;

    if (av < bv) {
      hi = mid;
    } else if (av == bv) {  // when av and bv are equal, need to move along a,
                            // so bring down the high search index along diag
      hi = mid;
    } else {  // av > bv
      lo = mid + 1;
    }
  }

  res.ai = diag <= a.size() ? diag - hi : a.size() - hi;
  res.bi = diag <= a.size() ? hi : diag - a.size() + hi;

  // {
  //     DiagonalSearchResult<a_view_t, b_view_t> res2 =
  //     diagonal_search2(a,b,diag); if (res.ai != res2.ai || res.bi != res2.bi)
  //     {
  //         printf("ai diag=%d expected=%d,%d actual=%d,%d\n",
  //         int(diag), int(res.ai), int(res2.ai), int(res.bi), int(res2.bi));
  //     }
  // }

  return res;
}

/*
 */
template <typename View>
KOKKOS_INLINE_FUNCTION DiagonalSearchResult<typename View::size_type, typename View::size_type> diagonal_search(
    const View &a, typename View::non_const_value_type totalWork, typename View::size_type diag) {
  using value_type = typename View::non_const_value_type;
  using size_type  = typename View::size_type;

  KokkosKernels::Impl::Iota<value_type, size_type> iota(totalWork);
  return diagonal_search(a, iota, diag);
}

}  // namespace Impl
}  // namespace KokkosGraph

#endif  // _KOKKOSGRAPH_MERGEPATH_IMPL_HPP
