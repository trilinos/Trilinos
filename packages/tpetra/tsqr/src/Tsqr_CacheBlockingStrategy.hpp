// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_CacheBlockingStrategy_hpp
#define __TSQR_CacheBlockingStrategy_hpp

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility> // std::pair

namespace TSQR {

  /// \class CacheBlockingStrategy
  /// \brief Tells CacheBlocker how to block up a tall skinny matrix.
  /// \author Mark Hoemmen
  ///
  /// This "strategy object" helps CacheBlocker decide how to block up
  /// a given tall skinny matrix by row into cache blocks.  It knows
  /// how to find the location (row index) and number of rows of any
  /// cache block in the matrix.  You can use this either for
  /// parallelization (e.g., partitioning the matrix among processors
  /// in a way that respects cache blocks) or for \c SequentialTsqr
  /// (whose factor() routine iterates top-down through cache blocks,
  /// and whose apply() and explicit_Q() routines iterate bottom-up
  /// through cache blocks).
  ///
  /// The cache blocking strategy is formulated in terms of \c
  /// SequentialTsqr.  All intranode parallel TSQR implementations use
  /// either SequentialTsqr or an algorithm like it, so this
  /// formulation is general enough for our needs.
  template<class LocalOrdinal, class Scalar>
  class CacheBlockingStrategy {
  public:
    /// \brief Constructor
    ///
    /// The cache blocking strategy asks for a cache size hint.  The
    /// appropriate cache level to use depends on the bandwidth of
    /// each cache level, whether it is shared among cores, and other
    /// hardware-specific features that are hard for us to model or
    /// measure.  In general, though, if each CPU core has its own L2
    /// cache, it would be appropriate to use that cache's size.
    ///
    /// The cache size need not be exact, but picking too small or too
    /// large a value will make TSQR slower.  In practice, TSQR is not
    /// sensitive to this value.  The sizeOfScalar parameter affects
    /// performance, not correctness (more or less -- it should never
    /// be zero, for example).  It's OK for it to be a slight
    /// overestimate.  Being much too big may affect performance in
    /// the same way as an excessively small cache size hint, and
    /// being much too small may affect performance in the same way as
    /// an excessively big cache size hint.
    ///
    /// \param cacheSizeHint [in] Cache size hint in bytes.  This is
    ///   used to pick the number of rows in a cache block.  If zero,
    ///   we guess a reasonable value.  This is a hint only, not a
    ///   command; the strategy may revise this, but it will not
    ///   change the revised value (that is, \c cache_size_hint() is a
    ///   constant for this instance's lifetime).
    ///
    /// \param sizeOfScalar [in] The number of bytes required to store
    ///   a Scalar value.  This is used to compute the dimensions of
    ///   cache blocks.  If sizeof(Scalar) correctly reports the size
    ///   of the representation of Scalar in memory, you can use the
    ///   default.  The default is correct for float, double, and any
    ///   of various fixed-length structs (like double-double and
    ///   quad-double).  It should also work for std::complex<T> where
    ///   T is anything in the previous sentence's list.  It does
    ///   <it>not</it> work for arbitrary-precision types whose
    ///   storage is dynamically allocated, even if the amount of
    ///   storage is a constant.  In the latter case, you should
    ///   specify a nondefault value.
    ///
    /// \note If Scalar is an arbitrary-precision type whose
    ///   representation length can change at runtime, you should
    ///   construct a new CacheBlockingStrategy object whenever the
    ///   representation length changes.
    CacheBlockingStrategy (const size_t cacheSizeHint = 0,
                           const size_t sizeOfScalar = sizeof(Scalar)) :
      size_of_scalar_ (sizeOfScalar),
      cache_size_hint_ (default_cache_size_hint (cacheSizeHint, sizeOfScalar))
    {}

    //! Copy constructor
    CacheBlockingStrategy (const CacheBlockingStrategy& rhs) :
      size_of_scalar_ (rhs.size_of_scalar_),
      cache_size_hint_ (rhs.cache_size_hint())
    {}

    //! Assignment operator
    CacheBlockingStrategy& operator= (const CacheBlockingStrategy& rhs) {
      size_of_scalar_ = rhs.size_of_scalar_;
      cache_size_hint_ = rhs.cache_size_hint();
      return *this;
    }

    /// \brief The cache size hint in bytes.
    ///
    /// This may not necessarily equal the suggested cache size (input
    /// to the constructor).  We treat that as a hint rather than a
    /// command.
    ///
    /// \note It may make sense to vary the cache size hint at run
    ///   time, for automatic performance tuning (trying to guess the
    ///   optimal value) or adaptivity to varying load.  However, the
    ///   cache blocking strategy object must not change the cache
    ///   size hint during the object's lifetime, since that would
    ///   prevent correct manipulation of matrices with contiguously
    ///   stored cache blocks.  This is because cache block parameters
    ///   depend on the cache size hint.  Thus, client code may assume
    ///   that this method always returns the same value for the
    ///   lifetime of the strategy object.
    size_t cache_size_hint () const { return cache_size_hint_; }

    /// \brief Size of a Scalar object in bytes.
    ///
    /// See the constructor documentation for an explanation of why
    /// this may not necessarily be sizeof(Scalar).  It should be in
    /// most cases, however.
    size_t size_of_scalar () const { return size_of_scalar_; }

    /// \brief Pointer offset for the cache block with the given index.
    ///
    /// The pointer offset depends on whether cache blocks are stored
    /// contiguously in the matrix.  If the cache block index is out
    /// of range, the returned result is undefined.
    ///
    /// \param index [in] Zero-based index of the cache block.
    /// \param nrows [in] Number of rows in the matrix.
    /// \param ncols [in] Number of columns in the matrix.
    /// \param nrows_cache_block [in] The value returned by
    ///   \c cache_block_num_rows().
    /// \param contiguous_cache_blocks [in] Whether the cache
    ///   blocks in the matrix are stored contiguously.
    LocalOrdinal
    cache_block_offset (const LocalOrdinal index,
                        const LocalOrdinal nrows,
                        const LocalOrdinal ncols,
                        const LocalOrdinal nrows_cache_block,
                        const bool contiguous_cache_blocks) const
    {
      // Suppress compiler warning for the unused argument.
      (void) nrows;

      const LocalOrdinal my_row_start = index * nrows_cache_block;
      if (contiguous_cache_blocks)
        return my_row_start * ncols;
      else // the common case
        return my_row_start;
    }

    /// \brief Leading dimension (a.k.a. stride) of the cache block.
    ///
    /// If cache blocks are stored contiguously, their leading
    /// dimension may vary.  Otherwise, their leading dimension is
    /// just that of the whole matrix.
    ///
    /// \param index [in] Zero-based index of the cache block.
    /// \param nrows [in] Number of rows in the matrix.
    /// \param ncols [in] Number of columns in the matrix.
    /// \param lda [in] Leading dimension (a.k.a. stride) of
    ///   the matrix.
    /// \param nrows_cache_block [in] The value returned by
    ///   \c cache_block_num_rows().
    /// \param contiguous_cache_blocks [in] Whether the cache
    ///   blocks in the matrix are stored contiguously.
    LocalOrdinal
    cache_block_stride (const LocalOrdinal index,
                        const LocalOrdinal nrows,
                        const LocalOrdinal ncols,
                        const LocalOrdinal lda,
                        const LocalOrdinal nrows_cache_block,
                        const bool contiguous_cache_blocks) const
    {
      if (contiguous_cache_blocks) {
        std::pair<LocalOrdinal, LocalOrdinal> result =
          cache_block (index, nrows, ncols, nrows_cache_block);
        return result.second; // Number of rows in the cache block
      }
      else {
        return lda;
      }
    }

    /// \brief Start and size of cache block number \c index.
    ///
    /// \param index [in] Zero-based index of the cache block.
    /// \param nrows [in] Number of rows in the whole matrix.
    /// \param ncols [in] Number of columns in the whole matrix.
    /// \param nrows_cache_block [in] The value returned by
    ///   \c cache_block_num_rows().
    ///
    /// \return If the input \c index is in range: The starting row
    ///   index (zero-based) of the cache block, and the number of
    ///   rows in the cache block.  If the input \c index is out of
    ///   range, then (nrows, 0).
    std::pair<LocalOrdinal, LocalOrdinal>
    cache_block (const LocalOrdinal index,
                 const LocalOrdinal nrows,
                 const LocalOrdinal ncols,
                 const LocalOrdinal nrows_cache_block) const
    {
      // See the comments in num_cache_blocks() for an explanation how
      // the number of cache blocks is computed, so that no cache
      // block has fewer than ncols rows.
      const LocalOrdinal quotient = nrows / nrows_cache_block;
      const LocalOrdinal remainder = nrows - nrows_cache_block * quotient;
      LocalOrdinal my_row_start, my_nrows;

      my_row_start = index * nrows_cache_block;
      if (quotient == 0) { // There is only one cache block.
        if (index == 0) {
          my_nrows = remainder;
        }
        else {
          my_nrows = 0; // Out-of-range block, therefore empty
        }
      }
      else if (remainder < ncols) { // There are quotient cache blocks.
        if (index < 0) {
          my_nrows = 0; // Out-of-range block, therefore empty
        }
        else if (index < quotient - 1) {
          my_nrows = nrows_cache_block;
        }
        else if (index == quotient - 1) {
          // The last cache block gets the leftover rows, so that no
          // cache block has fewer than ncols rows.
          my_nrows = nrows_cache_block + remainder;
        }
        else {
          my_nrows = 0; // Out-of-range block, therefore empty
        }
      }
      else { // There are quotient+1 cache blocks.
        if (index < 0) {
          my_nrows = 0; // Out-of-range block, therefore empty
        }
        else if (index < quotient) {
          my_nrows = nrows_cache_block;
        }
        else if (index == quotient) {
          // The last cache block has the leftover rows, which are
          // >= ncols and < nrows_cache_block.
          my_nrows = remainder;
        }
        else {
          my_nrows = 0; // Out-of-range block, therefore empty
        }
      }
      return std::make_pair (my_row_start, my_nrows);
    }

    /// \brief Complete description of the cache block.
    ///
    /// "Complete" means that it includes the location as well as the
    /// layout of the cache block.  This lets you construct a view of
    /// the cache block right away, given a view of the whole matrix.
    ///
    /// \param index [in] Zero-based index of the cache block.
    /// \param nrows [in] Number of rows in the whole matrix.
    /// \param ncols [in] Number of columns in the whole matrix.
    /// \param lda [in] Leading dimension (a.k.a. stride) of
    ///   the whole matrix.
    /// \param nrows_cache_block [in] The value returned by
    ///   \c cache_block_num_rows().
    /// \param contiguous_cache_blocks [in] Whether the cache
    ///   blocks in the matrix are stored contiguously.
    ///
    /// \return Four LocalOrdinals: The starting row index, the number
    ///   of rows in the cache block, the pointer offset of the cache
    ///   block, and leading dimension of the cache block.
    ///
    /// \note This method has an \f$O(1)\f$ cost, so that
    ///   parallelization by calling this method repeatedly for a
    ///   sequence of cache block indices is not expensive.
    std::vector<LocalOrdinal>
    cache_block_details (const LocalOrdinal index,
                         const LocalOrdinal nrows,
                         const LocalOrdinal ncols,
                         const LocalOrdinal lda,
                         const LocalOrdinal nrows_cache_block,
                         const bool contiguous_cache_blocks) const
    {
      const std::pair<LocalOrdinal, LocalOrdinal> result =
        cache_block (index, nrows, ncols, nrows_cache_block);
      const LocalOrdinal my_row_start = result.first;
      const LocalOrdinal my_nrows = result.second;

      const LocalOrdinal offset =
        contiguous_cache_blocks ? my_row_start * ncols : my_row_start;
      const LocalOrdinal stride =
        contiguous_cache_blocks ? my_nrows : lda;

      std::vector<LocalOrdinal> retval (4);
      retval[0] = my_row_start;
      retval[1] = my_nrows;
      retval[2] = offset;
      retval[3] = stride;
      return retval;
    }

    /// \brief Total number of cache blocks in the matrix.
    ///
    /// The input matrix is nrows by ncols.  The suggested number of
    /// rows per cache block is nrows_cache_block, but some cache
    /// blocks may have more or less rows.  However, no cache block
    /// may have fewer than ncols rows.
    ///
    /// \param nrows [in] Number of rows in the matrix.
    /// \param ncols [in] Number of columns in the matrix.
    /// \param nrows_cache_block [in] The value returned by
    ///   \c cache_block_num_rows().
    ///
    /// \return Total number of cache blocks in the matrix.
    ///
    LocalOrdinal
    num_cache_blocks (const LocalOrdinal nrows,
                      const LocalOrdinal ncols,
                      const LocalOrdinal nrows_cache_block) const
    {
      const LocalOrdinal quotient = nrows / nrows_cache_block;
      const LocalOrdinal remainder = nrows - nrows_cache_block * quotient;

      if (quotient == 0)
        // If nrows < nrows_cache_block, then there is only one cache
        // block, which gets all the rows.
        return static_cast<LocalOrdinal>(1);
      else if (remainder < ncols)
        // Don't let the last cache block have fewer than ncols rows.
        // If it would, merge it with the cache block above it.
        return quotient;
      else
        // The last cache block has the leftover rows, which are >=
        // ncols and < nrows_cache_block.
        return quotient + 1;
    }

    /// \brief Number of rows in the top cache block.
    ///
    /// If we partition the nrows by ncols matrix A into [A_top;
    /// A_bot] with A_top being a cache block and A_bot being the rest
    /// of the matrix, return the number of rows that A_top should
    /// have.
    ///
    /// \param nrows [in] "Current" number of rows in the matrix.  We
    ///   write "current" because this method is meant to be called
    ///   recursively over the "rest" of the matrix, until the "rest"
    ///   of the matrix has no more rows (is "empty").
    ///
    /// \param ncols [in] Number of columns in the matrix.
    ///
    /// \param nrows_cache_block [in] The value returned by
    ///   \c cache_block_num_rows().
    ///
    /// \return # of rows in top cache block A_top
    LocalOrdinal
    top_block_split_nrows (const LocalOrdinal nrows,
                           const LocalOrdinal ncols,
                           const LocalOrdinal nrows_cache_block) const
    {
      // We want to partition the nrows by ncols matrix A into [A_top;
      // A_bot], where A_top has nrows_cache_block rows.  However, we
      // don't want A_bot to have less than ncols rows.  If it would,
      // then we partition A so that A_top has nrows rows and A_bot is
      // empty.
      if (nrows < nrows_cache_block + ncols)
        return nrows;
      else
        // Don't ask for a bigger cache block than there are rows in
        // the matrix left to process.
        return std::min (nrows_cache_block, nrows);
    }

    /// \brief Number of rows in the bottom cache block.
    ///
    /// If we partition the nrows by ncols matrix A into [A_top;
    /// A_bot] with A_bot being a cache block and A_top being the rest
    /// of the matrix, return the number of rows that A_bot should
    /// have.
    ///
    /// \param nrows [in] "Current" number of rows in the matrix.  We
    ///   write "current" because this method is meant to be called
    ///   recursively over the "rest" of the matrix, until the "rest"
    ///   of the matrix has no more rows (is "empty").
    /// \param ncols [in] Number of columns in the matrix.
    /// \param nrows_cache_block [in] The value returned by
    ///   \c cache_block_num_rows().
    ///
    /// \return # of rows in top cache block A_bot
    LocalOrdinal
    bottom_block_split_nrows (const LocalOrdinal nrows,
                              const LocalOrdinal ncols,
                              const LocalOrdinal nrows_cache_block) const
    {
      // We split off the bottom block using the same splitting as if
      // we had split off as many top blocks of nrows_cache_block rows
      // as permissible.  The last block may have fewer than
      // nrows_cache_block rows, but it may not have fewer than ncols
      // rows (since we don't want any cache block to have fewer rows
      // than columns).
      const LocalOrdinal quotient = nrows / nrows_cache_block;
      const LocalOrdinal remainder = nrows - quotient * nrows_cache_block;

      LocalOrdinal nrows_bottom;
      if (quotient == 0)
        nrows_bottom = remainder;
      else if (remainder < ncols)
        nrows_bottom = nrows_cache_block + remainder;
      else if (remainder >= ncols)
        nrows_bottom = remainder;
      else
        throw std::logic_error("Should never get here!");
      return nrows_bottom;
    }

    /// \brief Default or revised cache size hint in bytes.
    ///
    /// If the input is zero, return a default cache size in bytes.
    /// Otherwise, revise the given suggestion based on the size of
    /// the Scalar type.  The result need not equal the suggested
    /// cache size, even if the latter is nonzero.  Call \c
    /// cache_size_hint() after calling this method, in order to get
    /// the actual cache size that the cache blocking strategy will
    /// use.
    ///
    /// \param suggested_cache_size [in] Suggested size of the cache
    ///   in bytes.  A hint, not a command.
    /// \param sizeOfScalar [in] Size of the Scalar type in bytes.
    ///
    /// \return Default or revised cache size hint in bytes.
    size_t
    default_cache_size_hint (const size_t suggested_cache_size,
                             const size_t sizeOfScalar) const
    {
      // This is a somewhat arbitrary minimum.  However, our TSQR
      // implementation was optimized for matrices with 20 or fewer
      // columns, and we expect matrices with 10 columns.  Thus, it's
      // reasonable to base the minimum on requiring that the cache
      // blocks for a matrix with 10 columns have no fewer rows than
      // columns.  We base the minimum on explicit_Q() with such a
      // matrix: a cache block of Q (10 x 10), a cache block of C (10
      // x 10), a TAU array (length 10), and the top block of C
      // (C_top) (10 x 10 in this case, and n x n in general for a
      // matrix with n columns).  In this bound, the cache blocks are
      // only square because of the requirement that they have no
      // fewer rows than columns; normally, cache blocks have many
      // more rows than columns.
      const size_t min_cache_size = sizeOfScalar * (3*10*10 + 10);

      // 64 KB is a reasonable guess for the L2 cache size.  If Scalar
      // is huge, min_cache_size above might be bigger, so we account
      // for that with a max.  The type cast is necessary so that the
      // compiler can decide which specialization of std::max to use
      // (the one for size_t, in this case, rather than the one for
      // int, which is the implicit type of the integer constant).
      const size_t default_cache_size =
        std::max (min_cache_size, static_cast<size_t> (65536));

      // If the suggested cache size is less than the minimum, ignore
      // the suggestion and pick the minimum.
      return (suggested_cache_size == 0) ?
        default_cache_size : std::max (suggested_cache_size, min_cache_size);
    }

    /// \brief "Typical" number of rows per cache block.
    ///
    /// This is a function of the cache block size and the number of
    /// columns in the matrix.  Not all cache blocks (in particular,
    /// the last one) will have this number of rows, but "most" will
    /// (hence "typical").  The returned value applies to \c
    /// SequentialTsqr::factor(), \c SequentialTsqr::apply(), and \c
    /// SequentialTsqr::explicit_Q().  In particular, we choose the
    /// number of rows per cache block so that when applying the
    /// implicitly stored Q factor (returned by \c factor()) to a
    /// matrix C with the same number of columns as Q was on input to
    /// \c factor(), then two cache blocks (one of Q and the other of
    /// C) will fit in cache.  This is the typical case when using
    /// TSQR to orthogonalize vectors.
    ///
    /// \param ncols [in] Number of columns in the matrix whose QR
    ///   factorization is to be computed using an intranode TSQR
    ///   implementation.
    LocalOrdinal
    cache_block_num_rows (const LocalOrdinal ncols) const
    {
      // Suppose the cache can hold W words (of size size_of_scalar_
      // bytes each).  We have to use the same number of rows per
      // cache block for both the factorization and applying the Q
      // factor.
      //
      // The factorization requires a working set of
      // ncols*(nrows_cache_block + ncols) + 2*ncols words:
      //
      // 1. One ncols by ncols R factor (not packed)
      // 2. One nrows_cache_block by ncols cache block
      // 3. tau array of length ncols
      // 4. workspace array of length ncols
      //
      // That means nrows_cache_block should be <= W/ncols - ncols - 2.
      //
      // Applying the Q factor to a matrix C with the same number of
      // columns as Q requires a working set of
      // 2*nrows_cache_block*ncols + ncols*ncols + 2*ncols
      //
      // 1. Cache block of Q: nrows_cache_block by ncols
      // 2. C_top block: ncols by ncols
      // 3. C_cur block: nrows_cache_block by ncols
      // 4. tau array of length ncols
      // 5. workspace array of length ncols
      //
      // That means nrows_cache_block should be <= (W/(2*N) - N/2 -
      // 1).  Obviously this is smaller than for the factorization, so
      // we use this formula to pick nrows_cache_block.  It should also
      // be at least ncols.

      const size_t W = cache_size_hint() / size_of_scalar_;

      // Compute everything in size_t first, and cast to LocalOrdinal
      // at the end.  This may avoid overflow if the cache is very
      // large and/or LocalOrdinal is very small (say a short int).
      //
      // Also, since size_t is unsigned, make sure that the
      // subtractions don't make it negative.  If it does, then either
      // ncols is too big or the cache is too small.

      const size_t term1 = W / (2*ncols);
      const size_t term2 = ncols / 2 + 1;

      if (term1 <= term2)
        {
          // The cache must be very small.  Just make the cache blocks
          // square.  That will be inefficient, but wil result in
          // correct behavior.
          return ncols;
        }
      else
        {
          // The compiler can easily prove that term1 - term2 > 0,
          // since we've gotten to this point.  Of course that's
          // assuming that C++ compilers are smart...
          const size_t nrows_cache_block =
            std::max (term1 - term2, static_cast<size_t>(ncols));

          // Make sure that nrows_cache_block fits in a LocalOrdinal
          // type.  We do so by casting the size_t to a LocalOrdinal
          // and then back into a size_t.  This should work in the
          // typical case of LocalOrdinal=int, and also whenever
          // LocalOrdinal's binary representation has no more bits
          // than that of size_t.
          const LocalOrdinal nrows_cache_block_as_lo =
            static_cast<LocalOrdinal> (nrows_cache_block);
          if (static_cast<size_t>(nrows_cache_block_as_lo) != nrows_cache_block)
            {
              std::ostringstream os;
              os << "Error:  While deciding on the number of rows in a cache "
                "block for sequential TSQR, the decided-upon number of rows "
                 << nrows_cache_block << " does not fit in a LocalOrdinal "
                "type, whose max value is "
                 << std::numeric_limits<LocalOrdinal>::max() << ".";
              throw std::range_error (os.str());
            }
          else
            return static_cast<LocalOrdinal> (nrows_cache_block);
        }
    }

  private:
    /// \brief Size in bytes required to store one Scalar object.
    ///
    /// This comes first, before \c cache_size_hint_, because
    /// computing a default value for the latter requires knowing the
    /// size of the Scalar type.
    size_t size_of_scalar_;

    /// \brief Cache size (hint) in bytes.
    ///
    /// This should only be set as the return value of \c
    /// default_cache_size_hint(), since that method revises the input
    /// for reasonableness (in particular so that it is not too
    /// small).
    size_t cache_size_hint_;

  };
} // namespace TSQR

#endif // __TSQR_CacheBlockingStrategy_hpp
