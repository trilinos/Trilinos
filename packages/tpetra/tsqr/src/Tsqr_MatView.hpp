//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER

#ifndef __TSQR_Tsqr_MatView_hpp
#define __TSQR_Tsqr_MatView_hpp

// Define for bounds checking and other safety features, undefine for speed.
// #define TSQR_MATVIEW_DEBUG 1

#ifdef TSQR_MATVIEW_DEBUG
#  include <limits>
#endif // TSQR_MATVIEW_DEBUG
#include <sstream>
#include <stdexcept>

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
    const ptrdiff_t A_lda (A.lda());
    const ptrdiff_t ncols (A.extent(1));
    const ptrdiff_t B_lda (B.lda());
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

#ifdef TSQR_MATVIEW_DEBUG
  template<class Ordinal, class Scalar>
  class MatViewVerify {
  public:
    static void
    verify (const Ordinal num_rows,
            const Ordinal num_cols,
            const Scalar* const A,
            const Ordinal leading_dim)
    {
      using std::endl;

      bool good = true;
      std::ostringstream os;
      if (! std::numeric_limits<Ordinal>::is_integer) {
        good = false;
        os << "Error: Ordinal type must be an integer.";
      }
      if (std::numeric_limits<Ordinal>::is_signed) {
        if (num_rows < 0) {
          good = false;
          os << "Error: num_rows (= " << num_rows << ") < 0.";
        }
        if (num_cols < 0) {
          good = false;
          os << "Error: num_cols (= " << num_cols << ") < 0.";
        }
        if (leading_dim < 0) {
          good = false;
          os << "Error: leading_dim (= " << leading_dim << ") < 0.";
        }
      }
      if (leading_dim < num_rows) {
        good = false;
        os << "Error: leading_dim (= " << leading_dim << ") < num_rows (= "
           << num_rows << ").";
      }
      if (! good) {
        throw std::invalid_argument (os.str ());
      }
    }
  };
#endif // TSQR_MATVIEW_DEBUG


  // Forward declaration
  template<class Ordinal, class Scalar>
  class ConstMatView;

  // Forward declaration
  template<class Ordinal, class Scalar>
  class Matrix;

  /// \class MatView
  ///
  /// A read-and-write nonowning view of a column-oriented matrix.
  template<class Ordinal, class Scalar>
  class MatView {
  public:
    using scalar_type = Scalar;
    using ordinal_type = Ordinal;
    using pointer = Scalar*;
    using reference = Scalar&;

    MatView () = default;

    MatView (const ordinal_type num_rows,
             const ordinal_type num_cols,
             pointer const A,
             const ordinal_type leading_dim) :
      nrows_(num_rows),
      ncols_(num_cols),
      lda_(leading_dim),
      A_(A)
    {
#ifdef TSQR_MATVIEW_DEBUG
      MatViewVerify<ordinal_type, scalar_type>::
        verify (num_rows, num_cols, A, leading_dim);
#endif // TSQR_MATVIEW_DEBUG
    }

    MatView (const MatView& view) = default;
    MatView& operator= (const MatView& view) = default;
    MatView (MatView&& view) = default;
    MatView& operator= (MatView&& view) = default;

    reference
    operator() (const ordinal_type i,
                const ordinal_type j) const
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits<ordinal_type>::is_signed) {
        if (i < 0 || i >= extent(0)) {
          throw std::invalid_argument("Row range invalid");
        }
        else if (j < 0 || j >= extent(1)) {
          throw std::invalid_argument("Column range invalid");
        }
      }
      else {
        if (i >= extent(0)) {
          throw std::invalid_argument("Row range invalid");
        }
        else if (j >= extent(1)) {
          throw std::invalid_argument("Column range invalid");
        }
      }
      if (A_ == nullptr) {
        throw std::logic_error("Attempt to reference NULL data");
      }
#endif // TSQR_MATVIEW_DEBUG
      return A_[i + j*lda()];
    }

    constexpr ordinal_type extent(const int r) const noexcept {
      return r == 0 ? nrows_ : (r == 1 ? ncols_ : ordinal_type(0));
    }

    constexpr ordinal_type stride(const int r) const noexcept {
      return r == 0 ? ordinal_type(1) : (r == 1 ? lda_ : ordinal_type(0));
    }

    constexpr ordinal_type lda() const noexcept {
      return stride(1);
    }

    /// \note The function is const, only because returning A_ doesn't
    /// change any members of *this.  Of course one may use the
    /// resulting pointer to fiddle with entries in the matrix, but
    /// that doesn't affect the MatView's properties.
    pointer data() const { return A_; }
    bool empty() const { return extent(0) == 0 || extent(1) == 0; }

    /// Return a "row block" (submatrix of consecutive rows in the
    /// inclusive range [firstRow,lastRow]).
    MatView row_block (const ordinal_type firstRow,
                       const ordinal_type lastRow)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits<ordinal_type>::is_signed) {
        if (firstRow < 0 || firstRow > lastRow || lastRow >= extent(0)) {
          throw std::invalid_argument ("Row range invalid");
        }
      }
      else {
        if (firstRow > lastRow || lastRow >= extent(0)) {
          throw std::invalid_argument ("Row range invalid");
        }
      }
#endif // TSQR_MATVIEW_DEBUG
      return MatView (lastRow - firstRow + 1, extent(1), data() + firstRow, lda());
    }

    /// Split off and return the top cache block of nrows_top rows.
    /// Modify *this to be the "rest" of the matrix.
    ///
    /// \note Only use this method to split off a single cache block.
    ///   It breaks if you try to use it otherwise.
    ///
    /// \param nrows_top [in] Number of rows in the top block (which
    ///   this method returns)
    ///
    /// \param b_contiguous_blocks [in] Whether or not the entries of
    ///   the top block are stored contiguously in *this.  The default
    ///   is no (false).
    ///
    /// \return The top block of nrows_top rows.  Data is a shallow
    ///   copy of the data in *this.
    MatView
    split_top (const ordinal_type nrows_top,
               const bool b_contiguous_blocks = false)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits<ordinal_type>::is_signed && nrows_top < 0) {
        std::ostringstream os;
        os << "nrows_top (= " << nrows_top << ") < 0";
        throw std::invalid_argument (os.str());
      }
      else if (nrows_top > extent(0)) {
        std::ostringstream os;
        os << "nrows_top (= " << nrows_top << ") > nrows (= " << extent(0) << ")";
        throw std::invalid_argument (os.str());
      }
#endif // TSQR_MATVIEW_DEBUG

      pointer const A_top_ptr = data();
      pointer A_rest_ptr;
      const ordinal_type nrows_rest = extent(0) - nrows_top;
      ordinal_type lda_top, lda_rest;
      if (b_contiguous_blocks) {
        lda_top = nrows_top;
        lda_rest = nrows_rest;
        A_rest_ptr = A_top_ptr + nrows_top * extent(1);
      }
      else {
        lda_top = lda();
        lda_rest = lda();
        A_rest_ptr = A_top_ptr + nrows_top;
      }
      MatView A_top (nrows_top, extent(1), data(), lda_top);
      A_ = A_rest_ptr;
      nrows_ = nrows_rest;
      lda_ = lda_rest;

      return A_top;
    }

    /// Split off and return the bottom block.  Modify *this to be the
    /// "rest" of the matrix.
    MatView
    split_bottom (const ordinal_type nrows_bottom,
                  const bool b_contiguous_blocks = false)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits<ordinal_type>::is_signed && nrows_bottom < 0) {
        throw std::invalid_argument ("nrows_bottom < 0");
      }
      if (nrows_bottom > extent(0)) {
        throw std::invalid_argument ("nrows_bottom > nrows");
      }
#endif // TSQR_MATVIEW_DEBUG

      pointer const A_rest_ptr = data();
      pointer A_bottom_ptr;
      const ordinal_type nrows_rest = extent(0) - nrows_bottom;
      ordinal_type lda_bottom, lda_rest;
      if (b_contiguous_blocks) {
        lda_bottom = nrows_bottom;
        lda_rest = extent(0) - nrows_bottom;
        A_bottom_ptr = A_rest_ptr + nrows_rest * extent(1);
      }
      else {
        lda_bottom = lda();
        lda_rest = lda();
        A_bottom_ptr = A_rest_ptr + nrows_rest;
      }
      MatView A_bottom (nrows_bottom, extent(1), A_bottom_ptr, lda_bottom);
      A_ = A_rest_ptr;
      nrows_ = nrows_rest;
      lda_ = lda_rest;

      return A_bottom;
    }

    bool operator== (const MatView& rhs) const {
      return extent(0) == rhs.extent(0) && extent(1) == rhs.extent(1) &&
        lda() == rhs.lda() && data() == rhs.data();
    }

    bool operator!= (const MatView& rhs) const {
      return extent(0) != rhs.extent(0) || extent(1) != rhs.extent(1) ||
        lda() != rhs.lda() || data() != rhs.data();
    }

  private:
    ordinal_type nrows_ = 0;
    ordinal_type ncols_ = 0;
    ordinal_type lda_ = 0;
    pointer A_ = nullptr;
  };

  /// \class ConstMatView
  ///
  /// A read-only view of a column-oriented matrix.
  template<class Ordinal, class Scalar>
  class ConstMatView {
  public:
    using scalar_type = Scalar;
    using ordinal_type = Ordinal;
    using pointer = const Scalar*;

    ConstMatView () = default;

    /// \note g++ with -Wall wants A_ to be initialized after lda_,
    /// otherwise it emits a compiler warning.
    ConstMatView (const ordinal_type num_rows,
                  const ordinal_type num_cols,
                  const scalar_type* const A,
                  const ordinal_type leading_dim) :
      nrows_(num_rows),
      ncols_(num_cols),
      lda_(leading_dim),
      A_(A)
    {
#ifdef TSQR_MATVIEW_DEBUG
      MatViewVerify<ordinal_type, scalar_type>::
        verify (num_rows, num_cols, A, leading_dim);
#endif // TSQR_MATVIEW_DEBUG
    }

    ConstMatView (const ConstMatView&) = default;
    ConstMatView& operator= (const ConstMatView&) = default;
    ConstMatView (ConstMatView&&) = default;
    ConstMatView& operator= (ConstMatView&&) = default;

    const scalar_type&
    operator() (const ordinal_type i, const ordinal_type j) const
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits<ordinal_type>::is_signed) {
        if (i < 0 || i >= extent(0)) {
          throw std::invalid_argument("Row range invalid");
        }
        else if (j < 0 || j >= extent(1)) {
          throw std::invalid_argument("Column range invalid");
        }
      }
      else {
        if (i >= extent(0)) {
          throw std::invalid_argument("Row range invalid");
        }
        else if (j >= extent(1)) {
          throw std::invalid_argument("Column range invalid");
        }
      }
      if (A_ == nullptr) {
        throw std::logic_error("Attempt to reference NULL data");
      }
#endif // TSQR_MATVIEW_DEBUG
      return A_[i + j*lda()];
    }

    constexpr ordinal_type extent(const int r) const noexcept {
      return r == 0 ? nrows_ : (r == 1 ? ncols_ : ordinal_type(0));
    }

    ordinal_type lda() const { return lda_; }

    pointer data() const { return A_; }

    bool empty() const { return extent(0) == 0 || extent(1) == 0; }

    /// Return a "row block" (submatrix of consecutive rows in the
    /// inclusive range [firstRow,lastRow]).
    ConstMatView
    rowBlock (const ordinal_type firstRow,
              const ordinal_type lastRow) const
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (firstRow < 0 || lastRow >= extent(0)) {
        throw std::invalid_argument ("Row range invalid");
      }
#endif // TSQR_MATVIEW_DEBUG
      return ConstMatView (lastRow - firstRow + 1, extent(1),
                           data() + firstRow, lda());
    }

    /// \brief Split off and return the top block.  Modify *this to be
    ///   the "rest" of the matrix.
    ///
    /// \note Only use this method to split off a single cache block.
    ///   It breaks if you try to use it otherwise.
    ///
    /// \param nrows_top [in] Number of rows in the top block (which
    ///   this method returns)
    ///
    /// \param b_contiguous_blocks [in] Whether or not the entries of
    ///   the top block are stored contiguously in *this.  The default
    ///   is no (false).
    ///
    /// \return The top block of nrows_top rows.  Data is a shallow
    ///   copy of the data in *this.
    ConstMatView split_top (const ordinal_type nrows_top,
                            const bool b_contiguous_blocks = false)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits<ordinal_type>::is_signed && nrows_top < 0) {
        throw std::invalid_argument ("nrows_top < 0");
      }
      if (nrows_top > extent(0)) {
        throw std::invalid_argument ("nrows_top > nrows");
      }
#endif // TSQR_MATVIEW_DEBUG

      pointer const A_top_ptr = data();
      pointer A_rest_ptr;
      const ordinal_type nrows_rest = extent(0) - nrows_top;
      ordinal_type lda_top, lda_rest;
      if (b_contiguous_blocks) {
        lda_top = nrows_top;
        lda_rest = nrows_rest;
        A_rest_ptr = A_top_ptr + nrows_top * extent(1);
      }
      else {
        lda_top = lda();
        lda_rest = lda();
        A_rest_ptr = A_top_ptr + nrows_top;
      }
      ConstMatView A_top (nrows_top, extent(1), data(), lda_top);
      A_ = A_rest_ptr;
      nrows_ = nrows_rest;
      lda_ = lda_rest;

      return A_top;
    }

    /// \brief Split off and return the bottom block.  Modify *this to
    ///   be the "rest" of the matrix.
    ConstMatView
    split_bottom (const ordinal_type nrows_bottom,
                  const bool b_contiguous_blocks = false)
    {
#ifdef TSQR_MATVIEW_DEBUG
      if (std::numeric_limits<ordinal_type>::is_signed && nrows_bottom < 0) {
        throw std::invalid_argument ("nrows_bottom < 0");
      }
      if (nrows_bottom > extent(0)) {
        throw std::invalid_argument ("nrows_bottom > nrows");
      }
#endif // TSQR_MATVIEW_DEBUG

      pointer const A_rest_ptr = data();
      pointer A_bottom_ptr;
      const ordinal_type nrows_rest = extent(0) - nrows_bottom;
      ordinal_type lda_bottom, lda_rest;
      if (b_contiguous_blocks) {
        lda_bottom = nrows_bottom;
        lda_rest = extent(0) - nrows_bottom;
        A_bottom_ptr = A_rest_ptr + nrows_rest * extent(1);
      }
      else {
        lda_bottom = lda();
        lda_rest = lda();
        A_bottom_ptr = A_rest_ptr + nrows_rest;
      }
      ConstMatView A_bottom (nrows_bottom, extent(1), A_bottom_ptr, lda_bottom);
      A_ = A_rest_ptr;
      nrows_ = nrows_rest;
      lda_ = lda_rest;

      return A_bottom;
    }

    bool operator== (const ConstMatView& rhs) const {
      return extent(0) == rhs.extent(0) && extent(1) == rhs.extent(1) &&
        lda() == rhs.lda() && data() == rhs.data();
    }

    bool operator!= (const ConstMatView& rhs) const {
      return extent(0) != rhs.extent(0) || extent(1) != rhs.extent(1) ||
        lda() != rhs.lda() || data() != rhs.data();
    }

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
    const ordinal_type stride = tgt.lda();
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
    if (tgt_nrows != ptrdiff_t (src.extent (0)) ||
        tgt_ncols != ptrdiff_t (src.extent (1))) {
      std::ostringstream os;
      os << "TSQR::deep_copy: dimensions of tgt (output matrix) and "
        "src (input matrix) are not compatible.  tgt is "
         << tgt.extent (0) << " x " << tgt.extent (1) << ", but src "
        "is " << src.extent (0) << " x " << src.extent (1) << ".";
      throw std::invalid_argument (os.str ());
    }
    for (ptrdiff_t j = 0; j < tgt_ncols; ++j) {
      auto* const tgt_j = &tgt(0,j);
      const auto* const src_j = &src(0,j);
      for (ptrdiff_t i = 0; i < tgt_nrows; ++i) {
        tgt_j[i] = src_j[i];
      }
    }
  }
} // namespace TSQR


#endif // __TSQR_Tsqr_MatView_hpp
