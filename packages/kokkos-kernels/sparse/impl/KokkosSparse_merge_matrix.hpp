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

#ifndef KOKKOSSPARSE_MERGEMATRIX_HPP
#define KOKKOSSPARSE_MERGEMATRIX_HPP

#include <type_traits>

#include "KokkosKernels_Iota.hpp"
#include "KokkosKernels_LowerBound.hpp"
#include "KokkosKernels_Predicates.hpp"
#include "KokkosKernels_SafeCompare.hpp"

/// \file KokkosSparse_merge_matrix.hpp

namespace KokkosSparse::Impl {

// a joint index into a and b
template <typename AIndex, typename BIndex>
struct MergeMatrixPosition {
  using a_index_type = AIndex;
  using b_index_type = BIndex;

  AIndex ai;
  BIndex bi;
};

/*! \class MergeMatrixDiagonal
    \brief a view into the entries of the Merge Matrix along a diagonal

   @tparam AView Type of the input view a, must be rank 1
   @tparam BViewLike Type of the view-like object b, must be Kokkos::View or
  KokkosKernels::Iota Example merge matrix M of two arrays A (vertical) and B
  (horizontal), as seen in Odeh, Green, Mwassi, Shmueli, Birk Merge Path -
  Parallel Merging Made Simple 2012 M[i,j] = 1 iff A[i] > B[j] operator(k)
  returns A[i] > B[j] at the kth entry of the diagonal

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
  Index into a diagonal from the bottom-left.
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
  static_assert(Kokkos::is_view_v<BViewLike> || KokkosKernels::Impl::is_iota_v<BViewLike>,
                "MergeMatrixDiagonal BViewLike must be Kokkos::View or "
                "KokkosKernels::Iota");
  static_assert(BViewLike::rank == 1, "MergeMatrixDiagonal BViewLike must be rank 1");

  using execution_space = typename AView::execution_space;

  /**
   * Define the types for index and value of each view
   */
  using a_index_type = typename AView::size_type;
  using b_index_type = typename BViewLike::size_type;
  using a_value_type = typename AView::non_const_value_type;
  using b_value_type = typename BViewLike::non_const_value_type;

  using position_type = MergeMatrixPosition<a_index_type, b_index_type>;

  // implement bare minimum parts of the view interface
  enum { rank = 1 };
  using non_const_value_type = bool;  ///< Merge matrix entries are 0 or 1.

  using size_type =
      typename std::conditional<sizeof(typename AView::size_type) >= sizeof(typename BViewLike::size_type),
                                typename AView::size_type,
                                typename BViewLike::size_type>::type;  ///< The larger of the two view types' size_types

  /** \brief Initializes the view a and view-like object b and the diagonal.
   */
  KOKKOS_INLINE_FUNCTION
  MergeMatrixDiagonal(const AView &a, const BViewLike &b, const size_type diagonal) : a_(a), b_(b), d_(diagonal) {}
  MergeMatrixDiagonal() = default;

  /**
   * Computes the position along a and b for a given diagonal di
   *
   * @param di Current diagonal
   * @return The MatrixPosition corresponding to the current diagonal
   */
  KOKKOS_INLINE_FUNCTION
  position_type position(const size_type &di) const noexcept {
    position_type pos;
    if (0 == d_) {
      pos.ai = 0;
      pos.bi = 0;
      return pos;
    } else {
      pos = diag_to_a_b(di);
      pos.ai += 1;
      return pos;
    }
  }

  /**
   * Compares a[i] > b[j] along the diagonal at entry di
   *
   * @param di Current diagonal
   * @return True if a[i] > b[j], false otherwise
   */
  KOKKOS_INLINE_FUNCTION
  bool operator()(const size_type di) const {
    position_type pos = diag_to_a_b(di);

    if (pos.ai >= typename position_type::a_index_type(a_.size())) {
      return true;  // on the +a side out of matrix bounds is 1
    } else if (pos.bi >= typename position_type::b_index_type(b_.size())) {
      return false;  // on the +b side out of matrix bounds is 0
    } else {
      return KokkosKernels::Impl::safe_gt(a_(pos.ai), b_(pos.bi));
    }
  }

  /**
   * Returns the length of the diagonal
   *
   * @return Length of the diagonal
   */
  KOKKOS_INLINE_FUNCTION
  size_type size() const noexcept {
    if (d_ <= size_type(a_.size()) && d_ <= size_type(b_.size())) {
      return d_;
    } else if (d_ > size_type(a_.size()) && d_ > size_type(b_.size())) {
      // TODO: this returns nonsense if d_ happens to be outside the merge
      // matrix
      return a_.size() + b_.size() - d_;
    } else {
      return KOKKOSKERNELS_MACRO_MIN(a_.size(), b_.size());
    }
  }

 private:
  /**
   * Translates an index along the diagonal to indices into a_ and b_
   *
   * @param di Current diagonal
   * @return The corresponding MatrixPosition with indices into a_ and b_
   */
  KOKKOS_INLINE_FUNCTION
  position_type diag_to_a_b(const size_type &di) const noexcept {
    position_type res;
    res.ai = d_ < size_type(a_.size()) ? (d_ - 1) - di : a_.size() - 1 - di;
    res.bi = d_ < size_type(a_.size()) ? di : d_ + di - a_.size();
    return res;
  }

  AView a_;      ///< The a view
  BViewLike b_;  ///< The b view
  size_type d_;  ///< diagonal
};

/*! \brief Return the first index on diagonal \code diag
           in the merge matrix of \code a and \code b that is not 1
This is effectively a lower-bound search on the merge matrix diagonal
where the predicate is "equals 1"
*/
template <typename AView, typename BViewLike>
KOKKOS_INLINE_FUNCTION typename MergeMatrixDiagonal<AView, BViewLike>::position_type diagonal_search(
    const AView &a, const BViewLike &b, typename MergeMatrixDiagonal<AView, BViewLike>::size_type diag) {
  // unmanaged view types for a and b
  using um_a_view = Kokkos::View<typename AView::value_type *, typename AView::device_type, Kokkos::MemoryUnmanaged>;
  using um_b_view =
      Kokkos::View<typename BViewLike::value_type *, typename BViewLike::device_type, Kokkos::MemoryUnmanaged>;

  um_a_view ua(a.data(), a.size());

  // if BViewLike is an Iota, pass it on directly to MMD,
  // otherwise, create an unmanaged view of B
  using b_type = typename std::conditional<KokkosKernels::Impl::is_iota<BViewLike>::value, BViewLike, um_b_view>::type;

  using MMD = MergeMatrixDiagonal<um_a_view, b_type>;
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
  return mmd.position(idx);
}

template <typename TeamMember, typename AView, typename BViewLike>
KOKKOS_INLINE_FUNCTION typename MergeMatrixDiagonal<AView, BViewLike>::position_type diagonal_search(
    const TeamMember &handle, const AView &a, const BViewLike &b,
    typename MergeMatrixDiagonal<AView, BViewLike>::size_type diag) {
  // unmanaged view types for a and b
  using um_a_view = Kokkos::View<typename AView::value_type *, typename AView::device_type, Kokkos::MemoryUnmanaged>;
  using um_b_view =
      Kokkos::View<typename BViewLike::value_type *, typename BViewLike::device_type, Kokkos::MemoryUnmanaged>;

  um_a_view ua(a.data(), a.size());

  // if BViewLike is an Iota, pass it on directly to MMD,
  // otherwise, create an unmanaged view of B
  using b_type = typename std::conditional<KokkosKernels::Impl::is_iota<BViewLike>::value, BViewLike, um_b_view>::type;

  using MMD = MergeMatrixDiagonal<um_a_view, b_type>;
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
  return mmd.position(idx);
}

/*! \brief

    \return A MergeMatrixDiagonal::position_type
 */
template <typename View>
KOKKOS_INLINE_FUNCTION auto diagonal_search(const View &a, typename View::non_const_value_type totalWork,
                                            typename View::size_type diag) {
  using value_type = typename View::non_const_value_type;
  using size_type  = typename View::size_type;

  KokkosKernels::Impl::Iota<value_type, size_type> iota(totalWork);
  return diagonal_search(a, iota, diag);
}

}  // namespace KokkosSparse::Impl

#endif  // KOKKOSSPARSE_MERGEMATRIX_HPP
