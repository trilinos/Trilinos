// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_MATVIEW_HPP
#define TSQR_MATVIEW_HPP

#include "Teuchos_TestForException.hpp"
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace TSQR {

  template<class Ordinal, class Scalar>
  class MatView;

  template<class LO, class SC, class SourceScalar>
  void
  deep_copy (const MatView<LO, SC>& tgt,
             const SourceScalar& src);

  template<class TargetOrdinal, class TargetScalar,
           class SourceOrdinal, class SourceScalar,
           template<class LO, class SC> class SourceMat>
  void
  deep_copy (const MatView<TargetOrdinal, TargetScalar>& tgt,
             const SourceMat<SourceOrdinal, SourceScalar>& src);

  template<class FirstMatrixViewType, class SecondMatrixViewType>
  bool
  matrix_equal (const FirstMatrixViewType& A,
                const SecondMatrixViewType& B)
  {
    if (A.extent(0) != B.extent(0) || A.extent(1) != B.extent(1)) {
      return false;
    }
    const ptrdiff_t nrows (A.extent(0));
    const ptrdiff_t A_lda (A.stride(1));
    const ptrdiff_t ncols (A.extent(1));
    const ptrdiff_t B_lda (B.stride(1));
    const auto* A_j = A.data();
    const auto* B_j = B.data();
    for (ptrdiff_t j = 0; j < ncols; ++j, A_j += A_lda, B_j += B_lda) {
      for (ptrdiff_t i = 0; i < nrows; ++i) {
        if (A_j[i] != B_j[i]) {
          return false;
        }
      }
    }
    return true;
  }

  // Forward declaration
  template<class Ordinal, class Scalar>
  class Matrix;

  /// \class MatView
  ///
  /// A read-and-write nonowning view of a column-oriented matrix.
  template<class Ordinal, class Scalar>
  class MatView {
  public:
    using non_const_value_type = typename std::remove_const<Scalar>::type;
    using const_value_type = const non_const_value_type;
    using ordinal_type = Ordinal;
    using pointer = Scalar*;
    using const_pointer = const Scalar*;
    using reference = Scalar&;
    using const_reference = const Scalar&;

    MatView () = default;

    MatView (const ordinal_type num_rows,
             const ordinal_type num_cols,
             pointer const A,
             const ordinal_type leading_dim) :
      nrows_(num_rows),
      ncols_(num_cols),
      lda_(leading_dim),
      A_(A)
    {}

    MatView (const MatView& view) = default;
    MatView& operator= (const MatView& view) = default;
    MatView (MatView&& view) = default;
    MatView& operator= (MatView&& view) = default;

    // Participates in overload resolution only if the type of
    // rhs.data() is assignable to A_.
    template<class InputScalarType>
    MatView (const MatView<Ordinal, InputScalarType>& rhs) :
      nrows_ (rhs.extent(0)),
      ncols_ (rhs.extent(1)),
      lda_ (rhs.stride(1)),
      A_ (rhs.data())
    {}

    constexpr ordinal_type extent(const int r) const noexcept {
      return r == 0 ? nrows_ : (r == 1 ? ncols_ : ordinal_type(0));
    }

    constexpr ordinal_type stride(const int r) const noexcept {
      return r == 0 ? ordinal_type(1) : (r == 1 ? lda_ : ordinal_type(0));
    }

    reference
    operator() (const ordinal_type i,
                const ordinal_type j) const
    {
      return A_[i + j * this->stride(1)];
    }

    pointer data() const { return A_; }

  private:
    ordinal_type nrows_ = 0;
    ordinal_type ncols_ = 0;
    ordinal_type lda_ = 0;
    pointer A_ = nullptr;
  };

  template<class LO, class SC, class SourceScalar>
  void
  deep_copy (const MatView<LO, SC>& tgt, const SourceScalar& src)
  {
    using ordinal_type = typename MatView<LO, SC>::ordinal_type;
    const ordinal_type num_rows = tgt.extent(0);
    const ordinal_type num_cols = tgt.extent(1);
    const ordinal_type stride = tgt.stride(1);
    auto* tgt_j = tgt.data();
    for (ordinal_type j = 0; j < num_cols; ++j, tgt_j += stride) {
      for (ordinal_type i = 0; i < num_rows; ++i) {
        tgt_j[i] = src;
      }
    }
  }

  template<class TargetOrdinal, class TargetScalar,
           class SourceOrdinal, class SourceScalar,
           template<class LO, class SC> class SourceMat>
  void
  deep_copy (const MatView<TargetOrdinal, TargetScalar>& tgt,
             const SourceMat<SourceOrdinal, SourceScalar>& src)
  {
    const ptrdiff_t tgt_nrows (tgt.extent (0));
    const ptrdiff_t tgt_ncols (tgt.extent (1));

    if (tgt_nrows == ptrdiff_t (src.extent (0)) ||
        tgt_ncols == ptrdiff_t (src.extent (1))) {
      for (ptrdiff_t j = 0; j < tgt_ncols; ++j) {
        auto* const tgt_j = &tgt(0,j);
        const auto* const src_j = &src(0,j);
        for (ptrdiff_t i = 0; i < tgt_nrows; ++i) {
          tgt_j[i] = src_j[i];
        }
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::invalid_argument, "TSQR::deep_copy: dimensions "
         "of tgt (output matrix) and src (input matrix) are not "
         "compatible.  tgt is " << tgt.extent (0) << " x " <<
         tgt.extent (1) << ", but src is " << src.extent (0) << " x "
         << src.extent (1) << ".");
    }
  }

  template<class MatViewType>
  std::pair<MatViewType, MatViewType>
  partition_2x1 (const MatViewType& A,
                 const typename MatViewType::ordinal_type nrows_top,
                 const bool b_contiguous_blocks = false)
  {
    using ordinal_type = typename MatViewType::ordinal_type;
    using pointer = typename MatViewType::pointer;

    const ordinal_type ncols = A.extent(1);
    pointer const A_top_ptr = A.data();
    const ordinal_type nrows_bot = A.extent(0) - nrows_top;

    pointer A_bot_ptr;
    ordinal_type lda_top, lda_bot;
    if (b_contiguous_blocks) {
      lda_top = nrows_top;
      lda_bot = nrows_bot;
      A_bot_ptr = A_top_ptr + nrows_top * A.extent(1);
    }
    else { // assume column major (LayoutLeft, in Kokkos terms)
      lda_top = A.stride(1);
      lda_bot = A.stride(1);
      A_bot_ptr = A_top_ptr + nrows_top;
    }

    MatViewType A_top (nrows_top, ncols, A_top_ptr, lda_top);
    MatViewType A_bot (nrows_bot, ncols, A_bot_ptr, lda_bot);
    return {A_top, A_bot};
  }

  template<class MatViewType>
  std::pair<MatViewType, MatViewType>
  partition_1x2 (const MatViewType& A,
                 const typename MatViewType::ordinal_type ncols_left)
  {
    using ordinal_type = typename MatViewType::ordinal_type;
    using pointer = typename MatViewType::pointer;

    const ordinal_type nrows = A.extent(0);
    const ordinal_type ncols = A.extent(1);
    const ordinal_type ncols_right = ncols - ncols_left;
    // assumes column major
    const auto right_offset = A.stride(1) * ncols_right;

    pointer A_top_ptr = A.data();
    pointer A_bot_ptr = A.data() + right_offset;

    MatViewType A_top (nrows, ncols_left, A_top_ptr, A.stride(1));
    MatViewType A_bot (nrows, ncols_right, A_bot_ptr, A.stride(1));
    return {A_top, A_bot};
  }

  /// \brief Split off and return the top block of nrows_top rows.
  ///   Modify A in place to be the "rest" of the matrix.
  ///
  /// \param A [in] On input: The whole matrix view.  On output: A
  ///   view of the "rest" of the matrix, that is, the part "below"
  ///   the returned matrix view.
  ///
  /// \param nrows_top [in] Number of rows in the top block (which
  ///   this method returns).
  ///
  /// \param contiguousCacheBlocks [in] Whether or not the entries of
  ///   the top block are stored contiguously in A.  The default is no
  ///   (false).
  ///
  /// \return A view of the top block of nrows_top rows.
  template<class LO, class SC>
  MatView<LO, SC>
  split_top (MatView<LO, SC>& A,
             const LO nrows_top,
             const bool contiguousCacheBlocks = false)
  {
    using pointer = typename MatView<LO, SC>::pointer;
    pointer A_top_ptr = A.data();
    pointer A_rest_ptr {};
    const LO nrows_rest = A.extent(0) - nrows_top;
    const LO ncols = A.extent(1);

    LO lda_top, lda_rest;
    if (contiguousCacheBlocks) {
      lda_top = nrows_top;
      lda_rest = nrows_rest;
      A_rest_ptr = A_top_ptr + nrows_top * ncols;
    }
    else {
      lda_top = A.stride(1);
      lda_rest = A.stride(1);
      A_rest_ptr = A_top_ptr + nrows_top;
    }
    MatView<LO, SC> A_top (nrows_top, ncols, A_top_ptr, lda_top);
    A = MatView<LO, SC> (nrows_rest, ncols, A_rest_ptr, lda_rest);
    return A_top;
  }

  /// \brief Split off and return the bottom block.  Modify A to be
  ///   the "rest" of the matrix.
  template<class LO, class SC>
  MatView<LO, SC>
  split_bottom (MatView<LO, SC>& A,
                const LO nrows_bottom,
                const bool contiguousCacheBlocks = false)
  {
    using pointer = typename MatView<LO, SC>::pointer;

    pointer A_rest_ptr = A.data();
    pointer A_bottom_ptr {};
    const LO nrows_rest = A.extent(0) - nrows_bottom;
    const LO ncols = A.extent(1);

    LO lda_bottom, lda_rest;
    if (contiguousCacheBlocks) {
      lda_bottom = nrows_bottom;
      lda_rest = A.extent(0) - nrows_bottom;
      A_bottom_ptr = A_rest_ptr + nrows_rest * ncols;
    }
    else {
      lda_bottom = A.stride(1);
      lda_rest = A.stride(1);
      A_bottom_ptr = A_rest_ptr + nrows_rest;
    }
    MatView<LO, SC> A_bottom (nrows_bottom, ncols, A_bottom_ptr, lda_bottom);
    A = MatView<LO, SC> (nrows_rest, ncols, A_rest_ptr, lda_rest);
    return A_bottom;
  }

  template<class LO, class SC>
  bool empty (const MatView<LO, SC>& A) {
    return A.extent(0) == 0 || A.extent(1) == 0;
  }

  template<class LO, class TargetScalar, class SourceScalar>
  void
  copy_upper_triangle (const MatView<LO, TargetScalar>& R_out,
                       const MatView<LO, SourceScalar>& R_in)
  {
    const LO nrows = R_out.extent (0);
    const LO ncols = R_out.extent (1);

    if (nrows >= ncols) {
      for (LO j = 0; j < ncols; ++j) {
        for (LO i = 0; i <= j; ++i) {
          R_out(i,j) = R_in(i,j);
        }
      }
    }
    else {
      auto R_out_lr = partition_1x2 (R_out, nrows);
      auto R_in_lr = partition_1x2 (R_in, nrows);
      copy_upper_triangle (R_out_lr.first, R_in_lr.first);
      deep_copy (R_out_lr.second, R_in_lr.second);
    }
  }

} // namespace TSQR

#endif // TSQR_MATVIEW_HPP
