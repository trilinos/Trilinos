// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_CacheBlocker_hpp
#define __TSQR_CacheBlocker_hpp

#include "Tsqr_CacheBlockingStrategy.hpp"
#include "Tsqr_MatView.hpp"
#include "Tsqr_Util.hpp"

#include <iterator>
#include <sstream>
#include <stdexcept>

namespace TSQR {

/// \class CacheBlocker
/// \brief Break a tall skinny matrix by rows into cache blocks.
/// \author Mark Hoemmen
///
/// A CacheBlocker uses a particular cache blocking strategy to
/// partition an nrows by ncols matrix by rows into cache blocks.
/// The entries in a cache block may be stored contiguously, or as
/// non-contiguous partitions of a matrix stored conventionally (in
/// column-major order).
///
/// The CacheBlocker blocks any matrix with the same number of rows
/// in the same way, regardless of the number of columns (the cache
/// blocking strategy's number of columns is set on construction).
/// This is useful for TSQR's apply() routine, which requires that
/// the output matrix C be blocked in the same way as the input
/// matrix Q (in which the Q factor is stored implicitly).
template <class Ordinal, class Scalar>
class CacheBlocker {
 private:
  using mat_view_type       = MatView<Ordinal, Scalar>;
  using const_mat_view_type = MatView<Ordinal, const Scalar>;

  void
  validate() {
    if (nrows_cache_block_ < ncols_) {
      std::ostringstream os;
      os << "The typical cache block size is too small.  Only "
         << nrows_cache_block_ << " rows fit, but every cache block needs "
                                  "at least as many rows as the number of columns "
         << ncols_
         << " in the matrix.";
      throw std::logic_error(os.str());
    }
  }

 public:
  /// \brief Constructor
  ///
  /// \param num_rows Number of rows in the matrix to block.
  /// \param num_cols Number of columns in the matrix to block.
  /// \param strategy Cache blocking strategy object (passed by copy).
  ///
  /// \note The CacheBlocker's number of columns may differ from the
  ///   number of columns associated with the cache blocking
  ///   strategy.  The strategy uses a fixed number of columns for
  ///   all matrices with the same number of rows, so that it blocks
  ///   all such matrices in the same way (at the same row indices).
  ///   This is useful for TSQR's apply() and explicit_Q() methods.
  CacheBlocker(const Ordinal num_rows,
               const Ordinal num_cols,
               const CacheBlockingStrategy<Ordinal, Scalar>& strategy)
    : nrows_(num_rows)
    , ncols_(num_cols)
    , strategy_(strategy)
    , nrows_cache_block_(strategy_.cache_block_num_rows(ncols_)) {
    validate();
  }

  //! Default constructor, so that CacheBlocker is DefaultConstructible.
  CacheBlocker() = default;

  //! Copy constructor
  CacheBlocker(const CacheBlocker& rhs)
    : nrows_(rhs.extent(0))
    , ncols_(rhs.extent(1))
    , strategy_(rhs.strategy_)
    , nrows_cache_block_(rhs.nrows_cache_block_) {}

  //! Assignment operator
  CacheBlocker& operator=(const CacheBlocker& rhs) {
    nrows_             = rhs.extent(0);
    ncols_             = rhs.extent(1);
    strategy_          = rhs.strategy_;
    nrows_cache_block_ = rhs.nrows_cache_block_;
    return *this;
  }

  //! Cache size hint (in bytes).
  size_t cache_size_hint() const { return strategy_.cache_size_hint(); }

  constexpr Ordinal extent(const int r) const noexcept {
    return r == 0 ? nrows_ : (r == 1 ? ncols_ : Ordinal(0));
  }

  /// \brief Split A in place into [A_top; A_rest].
  ///
  /// Return the topmost cache block A_top of A, and modify A in
  /// place to be the "rest" of the matrix A_rest.
  ///
  /// \param A [in/out] On input: view of the matrix to split.
  ///   On output: the "rest" of the matrix.  If there is only
  ///   one cache block, A_top contains all of the matrix and
  ///   A is empty on output.
  /// \param contiguous_cache_blocks [in] Whether cache blocks in
  ///   the matrix A are stored contiguously (default is false).
  ///
  /// \return View of the topmost cache block A_top.
  ///
  /// \note The number of rows in A_top depends on the number of
  ///   columns with which this CacheBlocker was set up (rather than
  ///   the number of columns in A, which may not be the same).  The
  ///   idea is to have the number and distribution of rows in the
  ///   cache blocks be the same as the original extent(0) by extent(1)
  ///   matrix with which this CacheBlocker was initialized.
  template <class MatrixViewType>
  MatrixViewType
  split_top_block(MatrixViewType& A,
                  const bool contiguous_cache_blocks) const {
    typedef typename MatrixViewType::ordinal_type ordinal_type;
    const ordinal_type nrows_top =
        strategy_.top_block_split_nrows(A.extent(0), extent(1),
                                        nrows_cache_block());
    return split_top(A, nrows_top, contiguous_cache_blocks);
  }

  /// \brief View of the topmost cache block of A.
  ///
  /// The matrix view A is copied so the view itself won't be modified.
  ///
  /// \param A [in] View of the matrix to block.
  /// \param contiguous_cache_blocks [in] Whether cache blocks in
  ///   the matrix A are stored contiguously (default is false).
  ///
  /// \return View of the topmost cache block of A.
  template <class MatrixViewType>
  MatrixViewType
  top_block(const MatrixViewType& A, const bool contiguous_cache_blocks) const {
    typedef typename MatrixViewType::ordinal_type ordinal_type;
    // Ignore the number of columns in A, since we want to block all
    // matrices using the same cache blocking strategy.
    const ordinal_type nrows_top =
        strategy_.top_block_split_nrows(A.extent(0), extent(1),
                                        nrows_cache_block());
    MatrixViewType A_copy(A);
    return split_top(A_copy, nrows_top, contiguous_cache_blocks);
  }

  /// \brief Split A in place into [A_rest; A_bot].
  ///
  /// Return the bottommost cache block A_bot of A, and modify A in
  /// place to be the "rest" of the matrix A_rest.
  ///
  /// \param A [in/out] On input: view of the matrix to split.  On
  ///   output: the "rest" of the matrix.  If there is only one
  ///   cache block, A_bot contains all of the matrix and A is empty
  ///   on output.
  /// \param contiguous_cache_blocks [in] Whether cache blocks in
  ///   the matrix A are stored contiguously (default is false).
  ///
  /// \return View of the bottommost cache block A_bot.
  ///
  template <class MatrixViewType>
  MatrixViewType
  split_bottom_block(MatrixViewType& A,
                     const bool contiguous_cache_blocks) const {
    typedef typename MatrixViewType::ordinal_type ordinal_type;
    // Ignore the number of columns in A, since we want to block all
    // matrices using the same cache blocking strategy.
    const ordinal_type nrows_bottom =
        strategy_.bottom_block_split_nrows(A.extent(0), extent(1),
                                           nrows_cache_block());
    // split_bottom() sets A to A_rest, and returns A_bot.
    return split_bottom(A, nrows_bottom, contiguous_cache_blocks);
  }

  /// \brief Fill the matrix A with zeros, respecting cache blocks.
  ///
  /// A specialization of this method for a particular
  /// MatrixViewType will only compile if MatrixViewType has a
  /// method "fill(const Scalar)" or "fill(const Scalar&)".  The
  /// intention is that the method be non-const and that it fill in
  /// the entries of the matrix with Scalar(0).
  ///
  /// \param A [in/out] View of the matrix to fill with zeros.
  ///
  /// \param contiguous_cache_blocks [in] Whether the cache blocks
  ///   in A are stored contiguously.
  ///
  template <class MatrixViewType>
  void
  fill_with_zeros(MatrixViewType A,
                  const bool contiguous_cache_blocks) const {
    // Note: if the cache blocks are stored contiguously, A.stride(1)
    // won't be the correct leading dimension of A, but it won't
    // matter: we only ever operate on A_cur here, and A_cur's
    // leading dimension is set correctly by split_top_block().
    while (!empty(A)) {
      // This call modifies the matrix view A, but that's OK since
      // we passed the input view by copy, not by reference.
      MatrixViewType A_cur = split_top_block(A, contiguous_cache_blocks);
      deep_copy(A_cur, Scalar{});
    }
  }

  /// \brief Fill the matrix A with zeros, respecting cache blocks.
  ///
  /// This version of the method takes a raw pointer and matrix
  /// dimensions, rather than a matrix view object.  If
  /// contiguous_cache_blocks==false, the matrix is stored either in
  /// column-major order with leading dimension lda; else, the
  /// matrix is stored in cache blocks, with each cache block's
  /// entries stored contiguously in column-major order.
  ///
  /// \param num_rows [in] Number of rows in the matrix A.
  /// \param num_cols [in] Number of columns in the matrix A.
  /// \param A [out] The matrix to fill with zeros.
  /// \param lda [in] Leading dimension (a.k.a. stride) of the
  ///   matrix A.
  /// \param contiguous_cache_blocks [in] Whether the cache blocks
  ///   in A are stored contiguously.
  void
  fill_with_zeros(const Ordinal num_rows,
                  const Ordinal num_cols,
                  Scalar A[],
                  const Ordinal lda,
                  const bool contiguous_cache_blocks) const {
    // We say "A_rest" because it points to the remaining part of
    // the matrix left to process; at the beginning, the "remaining"
    // part is the whole matrix, but that will change as the
    // algorithm progresses.
    //
    // Note: if the cache blocks are stored contiguously, lda won't
    // be the correct leading dimension of A, but it won't matter:
    // we only ever operate on A_cur here, and A_cur's leading
    // dimension is set correctly by split_top_block.
    mat_view_type A_rest(num_rows, num_cols, A, lda);

    while (!empty(A_rest)) {
      // This call modifies A_rest.
      mat_view_type A_cur = split_top_block(A_rest, contiguous_cache_blocks);
      deep_copy(A_cur, Scalar{});
    }
  }

  /// \brief Cache-block the given A_in matrix into A_out.
  ///
  /// Given an nrows by ncols (with nrows >= ncols) matrix A_in,
  /// stored in column-major order with leading dimension lda_in (>=
  /// nrows), copy it into A_out in a cache-blocked row block
  /// format.  Each cache block is a matrix in column-major order,
  /// and the elements of a cache block are stored consecutively in
  /// A_out.  The number of rows in each cache block depends on the
  /// cache-blocking strategy that this CacheBlocker uses.
  ///
  /// \param num_rows [in] Total number of rows in the matrices A_in and A_out
  /// \param num_cols [in] Number of columns in the matrices A_in and A_out
  /// \param A_out [out] nrows*ncols contiguous storage into which to write
  ///   the cache-blocked output matrix.
  /// \param A_in [in] nrows by ncols matrix, stored in column-major
  ///   order with leading dimension lda_in >= nrows
  /// \param lda_in [in] Leading dimension of the matrix A_in
  void
  cache_block(const Ordinal num_rows,
              const Ordinal num_cols,
              Scalar A_out[],
              const Scalar A_in[],
              const Ordinal lda_in) const {
    // We say "*_rest" because it points to the remaining part of
    // the matrix left to cache block; at the beginning, the
    // "remaining" part is the whole matrix, but that will change as
    // the algorithm progresses.
    const_mat_view_type A_in_rest(num_rows, num_cols, A_in, lda_in);
    // Leading dimension doesn't matter since A_out will be cache blocked.
    mat_view_type A_out_rest(num_rows, num_cols, A_out, lda_in);

    while (!empty(A_in_rest)) {
      if (empty(A_out_rest)) {
        throw std::logic_error("A_out_rest is empty, but A_in_rest is not");
      }
      // This call modifies A_in_rest.
      const_mat_view_type A_in_cur = split_top_block(A_in_rest, false);
      // This call modifies A_out_rest.
      mat_view_type A_out_cur = split_top_block(A_out_rest, true);
      deep_copy(A_out_cur, A_in_cur);
    }
  }

  //! "Un"-cache-block the given A_in matrix into A_out.
  void
  un_cache_block(const Ordinal num_rows,
                 const Ordinal num_cols,
                 Scalar A_out[],
                 const Ordinal lda_out,
                 const Scalar A_in[]) const {
    // We say "*_rest" because it points to the remaining part of
    // the matrix left to cache block; at the beginning, the
    // "remaining" part is the whole matrix, but that will change as
    // the algorithm progresses.
    //
    // Leading dimension doesn't matter since A_in is cache blocked.
    const_mat_view_type A_in_rest(num_rows, num_cols, A_in, lda_out);
    mat_view_type A_out_rest(num_rows, num_cols, A_out, lda_out);

    while (!empty(A_in_rest)) {
      if (empty(A_out_rest)) {
        throw std::logic_error("A_out_rest is empty, but A_in_rest is not");
      }
      // This call modifies A_in_rest.
      const_mat_view_type A_in_cur = split_top_block(A_in_rest, true);
      // This call modifies A_out_rest.
      mat_view_type A_out_cur = split_top_block(A_out_rest, false);
      deep_copy(A_out_cur, A_in_cur);
    }
  }

  /// \brief Return the cache block with index \c cache_block_index.
  ///
  /// \param A [in] The original matrix.
  /// \param cache_block_index [in] Zero-based index of the cache block.
  ///   If the index is out of bounds, silently return an empty matrix
  ///   view.
  /// \param contiguous_cache_blocks [in] Whether cache blocks are
  ///   stored contiguously.
  ///
  /// \return Cache block of A with the given index, or an empty
  ///   matrix view if the index is out of bounds.
  ///
  /// \note This method is templated on MatrixViewType, so that it
  ///   works with any matrix view type.
  template <class MatrixViewType>
  MatrixViewType
  get_cache_block(MatrixViewType A,
                  const typename MatrixViewType::ordinal_type cache_block_index,
                  const bool contiguous_cache_blocks) const {
    typedef typename MatrixViewType::ordinal_type ordinal_type;

    // Total number of cache blocks.
    const ordinal_type num_cache_blocks =
        strategy_.num_cache_blocks(A.extent(0), A.extent(1), nrows_cache_block());

    if (cache_block_index >= num_cache_blocks) {
      return MatrixViewType{};  // empty
    }
    // result[0] = starting row index of the cache block
    // result[1] = number of rows in the cache block
    // result[2] = pointer offset (A.data() + result[2])
    // result[3] = leading dimension (a.k.a. stride) of the cache block
    std::vector<Ordinal> result =
        strategy_.cache_block_details(cache_block_index, A.extent(0),
                                      A.extent(1), A.stride(1),
                                      nrows_cache_block(),
                                      contiguous_cache_blocks);
    if (result[1] == 0) {
      return MatrixViewType{};
    }

    // We expect that ordinal_type is signed, so adding signed
    // (ordinal_type) to unsigned (pointer) may raise compiler
    // warnings.
    return MatrixViewType(result[1], A.extent(1),
                          A.data() + size_t(result[2]),
                          result[3]);
  }

 private:
  //! Number of rows in the matrix to block.
  Ordinal nrows_ = 0;

  //! Number of columns in the matrix to block.
  Ordinal ncols_ = 0;

  //! Strategy used to break the matrix into cache blocks.
  CacheBlockingStrategy<Ordinal, Scalar> strategy_;

  /// \brief Number of rows in a "typical" cache block.
  ///
  /// We could instead use the strategy object to recompute this
  /// quantity each time, but we choose to cache the computed value
  /// here.  For an explanation of "typical," see the documentation
  /// of \c nrows_cache_block().
  Ordinal nrows_cache_block_ = 0;

  /// \brief Number of rows in a "typical" cache block.
  ///
  /// For an explanation of "typical," see the documentation of
  /// CacheBlockingStrategy.  In brief, some cache blocks may have
  /// more rows (up to but not including nrows_cache_block() +
  /// extent(1) rows), and some may have less (but no less than
  /// extent(1) rows).
  size_t nrows_cache_block() const { return nrows_cache_block_; }
};

/// \class CacheBlockRangeIterator
/// \brief Bidirectional iterator over a contiguous range of cache blocks.
/// \author Mark Hoemmen
///
/// "Contiguous range of cache blocks" means that the indices of the
/// cache blocks, as interpreted by the CacheBlocker object, are
/// contiguous.
template <class MatrixViewType>
class CacheBlockRangeIterator : public std::iterator<std::forward_iterator_tag, MatrixViewType> {
 public:
  using view_type    = MatrixViewType;
  using ordinal_type = typename MatrixViewType::ordinal_type;
  using scalar_type  = typename MatrixViewType::non_const_value_type;

  /// \brief Default constructor.
  ///
  /// \note To implementers: We only implement a default constructor
  ///   because all iterators (e.g., TrivialIterator) must be
  ///   DefaultConstructible.
  CacheBlockRangeIterator() = default;

  /// \brief Standard constructor.
  ///
  /// \param A [in] View of the matrix over whose cache block(s) to
  ///   iterate.
  /// \param strategy [in] Cache blocking strategy for a matrix with
  ///   the same number of rows as the matrix A.
  /// \param currentIndex [in] The iterator's current cache block index.
  /// \param reverse [in] Whether to iterate over the cache blocks
  ///   in reverse order of their indices.
  /// \param contiguousCacheBlocks [in] Whether cache blocks in the
  ///   matrix A are stored contiguously.
  CacheBlockRangeIterator(const MatrixViewType& A,
                          const CacheBlockingStrategy<ordinal_type, scalar_type>& strategy,
                          const ordinal_type currentIndex,
                          const bool reverse,
                          const bool contiguousCacheBlocks)
    : A_(A)
    , blocker_(A_.extent(0), A_.extent(1), strategy)
    , curInd_(currentIndex)
    , reverse_(reverse)
    , contiguousCacheBlocks_(contiguousCacheBlocks) {}

  //! Copy constructor.
  CacheBlockRangeIterator(const CacheBlockRangeIterator& rhs)
    : A_(rhs.A_)
    , blocker_(rhs.blocker_)
    , curInd_(rhs.curInd_)
    , reverse_(rhs.reverse_)
    , contiguousCacheBlocks_(rhs.contiguousCacheBlocks_) {}

  //! Assignment operator.
  CacheBlockRangeIterator& operator=(const CacheBlockRangeIterator& rhs) {
    A_                     = rhs.A_;
    blocker_               = rhs.blocker_;
    curInd_                = rhs.curInd_;
    reverse_               = rhs.reverse_;
    contiguousCacheBlocks_ = rhs.contiguousCacheBlocks_;
    return *this;
  }

  //! Prefix increment operator.
  CacheBlockRangeIterator& operator++() {
    if (reverse_)
      --curInd_;
    else
      ++curInd_;
    return *this;
  }

  /// \brief Postfix increment operator.
  ///
  /// This may be less efficient than prefix operator++, since the
  /// postfix operator has to make a copy of the iterator before
  /// modifying it.
  CacheBlockRangeIterator operator++(int) {
    CacheBlockRangeIterator retval(*this);
    operator++();
    return retval;
  }

  /// \brief Equality operator.
  ///
  /// Equality of cache block range iterators only tests the cache
  /// block index, not reverse-ness.  This means we can compare a
  /// reverse-direction iterator with a forward-direction iterator,
  /// and vice versa.
  bool operator==(const CacheBlockRangeIterator& rhs) {
    // Not correct, but fast.  Should return false for different A_
    // or different blocker_.
    return curInd_ == rhs.curInd_;
  }

  //! Inequality operator.
  bool operator!=(const CacheBlockRangeIterator& rhs) {
    // Not correct, but fast.  Should return false for different A_
    // or different blocker_.
    return curInd_ != rhs.curInd_;
  }

  /// \brief A view of the current cache block.
  ///
  /// If the current cache block index is invalid, this returns an
  /// empty cache block (that is, calling empty() on the returned
  /// view returns true).
  MatrixViewType operator*() const {
    return blocker_.get_cache_block(A_, curInd_, contiguousCacheBlocks_);
  }

 private:
  MatrixViewType A_;
  CacheBlocker<ordinal_type, scalar_type> blocker_;
  ordinal_type curInd_        = 0;
  bool reverse_               = false;
  bool contiguousCacheBlocks_ = false;
};

/// \class CacheBlockRange
/// \brief Collection of cache blocks with a contiguous range of indices.
/// \author Mark Hoemmen
///
/// We mean "collection" in the C++ sense: you can iterate over the
/// elements using iterators.  The iterators are valid only when the
/// CacheBlockRange is in scope, just like the iterators of
/// std::vector.
///
/// CacheBlockRange is useful for \c KokkosNodeTsqr, in particular
/// for \c FactorFirstPass and \c ApplyFirstPass.  Sequential TSQR's
/// factorization is forward iteration over the collection, and
/// applying the Q factor or computing the explicit Q factor is
/// iteration in the reverse direction (decreasing cache block
/// index).
///
/// This class is templated so that it works with any matrix view.
template <class MatrixViewType>
class CacheBlockRange {
 public:
  using view_type    = MatrixViewType;
  using ordinal_type = typename MatrixViewType::ordinal_type;
  using scalar_type  = typename MatrixViewType::non_const_value_type;

  //! Type of an iterator over the range of cache blocks.
  using iterator = CacheBlockRangeIterator<MatrixViewType>;

  /// \brief Constructor
  ///
  /// \param A [in] View of the matrix to factor.
  /// \param strategy [in] Cache blocking strategy to use (copied
  ///   on input).
  /// \param startIndex [in] Starting index of the cache block
  ///   sequence.
  /// \param endIndex [in] Ending index (exclusive) of the cache
  ///   block sequence.  Precondition: startIndex <= endIndex.  If
  ///   startIndex == endIndex, the sequence is empty.
  /// \param contiguousCacheBlocks [in] Whether cache blocks in the
  ///   matrix A are stored contiguously.
  CacheBlockRange(MatrixViewType A,
                  const CacheBlockingStrategy<ordinal_type, scalar_type>& strategy,
                  const ordinal_type startIndex,
                  const ordinal_type endIndex,
                  const bool contiguousCacheBlocks)
    : A_(A)
    , startIndex_(startIndex)
    , endIndex_(endIndex)
    , strategy_(strategy)
    , contiguousCacheBlocks_(contiguousCacheBlocks) {}

  bool empty() const {
    return startIndex_ >= endIndex_;
  }

  iterator begin() const {
    return iterator(A_, strategy_, startIndex_, false, contiguousCacheBlocks_);
  }

  iterator end() const {
    return iterator(A_, strategy_, endIndex_, false, contiguousCacheBlocks_);
  }

  iterator rbegin() const {
    return iterator(A_, strategy_, endIndex_ - 1, true, contiguousCacheBlocks_);
  }

  iterator rend() const {
    // Think about it: rbegin() == rend() means that rbegin() is invalid
    // and shouldn't be dereferenced.  rend() should never be dereferenced.
    return iterator(A_, strategy_, startIndex_ - 1, true, contiguousCacheBlocks_);
  }

 private:
  //! View of the matrix.
  MatrixViewType A_;

  /// \brief Starting index of the range of cache blocks.
  ///
  /// We always have startIndex_ <= endIndex_.  Reverse-order
  /// iteration is indicated by the iterator's reverse_ member
  /// datum.
  ordinal_type startIndex_;

  /// \brief Ending index (exclusive) of the range of cache blocks.
  ///
  /// See the documentation of startIndex_ for its invariant.
  ordinal_type endIndex_;

  //! Cache blocking strategy for a matrix with the same number of rows as A_.
  CacheBlockingStrategy<ordinal_type, scalar_type> strategy_;

  //! Whether the cache blocks of the matrix A_ are stored contiguously.
  bool contiguousCacheBlocks_;
};

}  // namespace TSQR

#endif  // __TSQR_CacheBlocker_hpp
