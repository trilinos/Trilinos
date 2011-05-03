//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_CacheBlocker_hpp
#define __TSQR_CacheBlocker_hpp

#include <Tsqr_CacheBlockingStrategy.hpp>
#include <Tsqr_MatView.hpp>
#include <Tsqr_Util.hpp>

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
  template<class Ordinal, class Scalar>
  class CacheBlocker {
  private:
    typedef MatView<Ordinal, Scalar> mat_view;
    typedef ConstMatView<Ordinal, Scalar> const_mat_view;

    void
    validate () 
    {
      if (nrows_cache_block_ < ncols_)
	{
	  std::ostringstream os;
	  os << "The typical cache block size is too small.  Only " 
	     << nrows_cache_block_ << " rows fit, but every cache block needs "
	    "at least as many rows as the number of columns " << ncols_ 
	     << " in the matrix.";
	  throw std::logic_error (os.str());
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
    CacheBlocker (const Ordinal num_rows,
		  const Ordinal num_cols,
		  const CacheBlockingStrategy<Ordinal, Scalar>& strategy) :
      nrows_ (num_rows), 
      ncols_ (num_cols), 
      strategy_ (strategy),
      nrows_cache_block_ (strategy_.cache_block_num_rows (ncols()))
    {
      validate ();
    }

    //! Copy constructor
    CacheBlocker (const CacheBlocker& rhs) :
      nrows_ (rhs.nrows()), 
      ncols_ (rhs.ncols()), 
      strategy_ (rhs.strategy_),
      nrows_cache_block_ (rhs.nrows_cache_block_)
    {}

    //! Assignment operator
    CacheBlocker& operator= (const CacheBlocker& rhs) {
      nrows_ = rhs.nrows();
      ncols_ = rhs.ncols();
      strategy_ = rhs.strategy_;
      nrows_cache_block_ = rhs.nrows_cache_block_;
      return *this;
    }

    //! Cache block size in bytes.
    size_t cache_block_size () const { return strategy_.cache_block_size(); }

    //! Number of rows in the matrix to block.
    Ordinal nrows () const { return nrows_; }

    //! Number of columns in the matrix to block.
    Ordinal ncols () const { return ncols_; }

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
    ///   cache blocks be the same as the original nrows() by ncols()
    ///   matrix with which this CacheBlocker was initialized.
    template< class MatrixViewType >
    MatrixViewType
    split_top_block (MatrixViewType& A, const bool contiguous_cache_blocks) const
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      const ordinal_type nrows_top = 
	strategy_.top_block_split_nrows (A.nrows(), ncols(), 
					 nrows_cache_block());
      // split_top() sets A to A_rest, and returns A_top.
      return A.split_top (nrows_top, contiguous_cache_blocks);
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
    template< class MatrixViewType >
    MatrixViewType
    top_block (const MatrixViewType& A, const bool contiguous_cache_blocks) const
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      // Ignore the number of columns in A, since we want to block all
      // matrices using the same cache blocking strategy.
      const ordinal_type nrows_top = 
	strategy_.top_block_split_nrows (A.nrows(), ncols(), 
					 nrows_cache_block());
      MatrixViewType A_copy (A);
      return A_copy.split_top (nrows_top, contiguous_cache_blocks);
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
    template< class MatrixViewType >
    MatrixViewType
    split_bottom_block (MatrixViewType& A, const bool contiguous_cache_blocks) const
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      // Ignore the number of columns in A, since we want to block all
      // matrices using the same cache blocking strategy.
      const ordinal_type nrows_bottom = 
	strategy_.bottom_block_split_nrows (A.nrows(), ncols(), 
					    nrows_cache_block());
      // split_bottom() sets A to A_rest, and returns A_bot.
      return A.split_bottom (nrows_bottom, contiguous_cache_blocks);
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
    template<class MatrixViewType>
    void
    fill_with_zeros (MatrixViewType A,
		     const bool contiguous_cache_blocks) const
    {
      // Note: if the cache blocks are stored contiguously, A.lda()
      // won't be the correct leading dimension of A, but it won't
      // matter: we only ever operate on A_cur here, and A_cur's
      // leading dimension is set correctly by split_top_block().
      while (! A.empty())
	{
	  // This call modifies the matrix view A, but that's OK since
	  // we passed the input view by copy, not by reference.
	  MatrixViewType A_cur = split_top_block (A, contiguous_cache_blocks);
	  A_cur.fill (Scalar(0));
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
    fill_with_zeros (const Ordinal num_rows,
		     const Ordinal num_cols,
		     Scalar A[],
		     const Ordinal lda, 
		     const bool contiguous_cache_blocks) const
    {
      // We say "A_rest" because it points to the remaining part of
      // the matrix left to process; at the beginning, the "remaining"
      // part is the whole matrix, but that will change as the
      // algorithm progresses.
      //
      // Note: if the cache blocks are stored contiguously, lda won't
      // be the correct leading dimension of A, but it won't matter:
      // we only ever operate on A_cur here, and A_cur's leading
      // dimension is set correctly by A_rest.split_top().
      mat_view A_rest (num_rows, num_cols, A, lda);

      while (! A_rest.empty())
	{
	  // This call modifies A_rest.
	  mat_view A_cur = split_top_block (A_rest, contiguous_cache_blocks);
	  A_cur.fill (Scalar(0));
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
    cache_block (const Ordinal num_rows,
		 const Ordinal num_cols,
		 Scalar A_out[],
		 const Scalar A_in[],
		 const Ordinal lda_in) const
    {
      // We say "*_rest" because it points to the remaining part of
      // the matrix left to cache block; at the beginning, the
      // "remaining" part is the whole matrix, but that will change as
      // the algorithm progresses.
      const_mat_view A_in_rest (num_rows, num_cols, A_in, lda_in);
      // Leading dimension doesn't matter since A_out will be cache blocked.
      mat_view A_out_rest (num_rows, num_cols, A_out, lda_in);

      while (! A_in_rest.empty())
	{
	  if (A_out_rest.empty())
	    throw std::logic_error("A_out_rest is empty, but A_in_rest is not");

	  // This call modifies A_in_rest.
	  const_mat_view A_in_cur = split_top_block (A_in_rest, false);

	  // This call modifies A_out_rest.
	  mat_view A_out_cur = split_top_block (A_out_rest, true);

	  copy_matrix (A_in_cur.nrows(), num_cols, A_out_cur.get(), 
		       A_out_cur.lda(), A_in_cur.get(), A_in_cur.lda());
	}
    }

    //! "Un"-cache-block the given A_in matrix into A_out.
    void
    un_cache_block (const Ordinal num_rows,
		    const Ordinal num_cols,
		    Scalar A_out[],
		    const Ordinal lda_out,		    
		    const Scalar A_in[]) const
    {
      // We say "*_rest" because it points to the remaining part of
      // the matrix left to cache block; at the beginning, the
      // "remaining" part is the whole matrix, but that will change as
      // the algorithm progresses.
      //
      // Leading dimension doesn't matter since A_in is cache blocked.
      const_mat_view A_in_rest (num_rows, num_cols, A_in, lda_out);
      mat_view A_out_rest (num_rows, num_cols, A_out, lda_out);

      while (! A_in_rest.empty())
	{
	  if (A_out_rest.empty())
	    throw std::logic_error("A_out_rest is empty, but A_in_rest is not");

	  // This call modifies A_in_rest.
	  const_mat_view A_in_cur = split_top_block (A_in_rest, true);

	  // This call modifies A_out_rest.
	  mat_view A_out_cur = split_top_block (A_out_rest, false);

	  copy_matrix (A_in_cur.nrows(), num_cols, A_out_cur.get(), 
		       A_out_cur.lda(), A_in_cur.get(), A_in_cur.lda());
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
    ///   works with both MatView and ConstMatView.
    template<class MatrixViewType>
    MatrixViewType 
    get_cache_block (MatrixViewType A,
		     const typename MatrixViewType::ordinal_type cache_block_index,
		     const bool contiguous_cache_blocks) const
    {
      using Teuchos::Tuple;
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      typedef typename MatrixViewType::scalar_type scalar_type;

      // Total number of cache blocks.
      const ordinal_type num_cache_blocks = 
	strategy_.num_cache_blocks (A.nrows(), A.ncols(), nrows_cache_block());

      if (cache_block_index >= num_cache_blocks)
	return MatrixViewType (0, 0, NULL, 0); // empty

      // result[0] = starting row index of the cache block
      // result[1] = number of rows in the cache block
      // result[2] = pointer offset (A.get() + result[2])
      // result[3] = leading dimension (a.k.a. stride) of the cache block
      Tuple<Ordinal, 4> result = 
	strategy_.cache_block_details (cache_block_index, A.nrows(), A.ncols(),
				       A.lda(), nrows_cache_block(), 
				       contiguous_cache_blocks);
      if (result[1] == 0)
	// For some reason, the cache block is empty.  	
	return MatrixViewType (0, 0, NULL, 0);

      // We expect that ordinal_type is signed, so adding signed
      // (ordinal_type) to unsigned (pointer) may raise compiler
      // warnings.
      return MatrixViewType (result[1], A.ncols(), 
			     A.get() + static_cast<size_t>(result[2]), 
			     result[3]);
    }

  private:
    //! Number of rows in the matrix to block.
    Ordinal nrows_;

    //! Number of columns in the matrix to block.
    Ordinal ncols_;
    
    //! Strategy used to break the matrix into cache blocks.
    CacheBlockingStrategy<Ordinal, Scalar> strategy_;

    /// \brief Number of rows in a "typical" cache block.
    ///
    /// We could instead use the strategy object to recompute this
    /// quantity each time, but we choose to cache the computed value
    /// here.  For an explanation of "typical," see the documentation
    /// of \c nrows_cache_block().
    Ordinal nrows_cache_block_;

    /// \brief Number of rows in a "typical" cache block.
    ///
    /// For an explanation of "typical," see the documentation of
    /// CacheBlockingStrategy.  In brief, some cache blocks may have
    /// more rows (up to but not including nrows_cache_block() +
    /// ncols() rows), and some may have less (but no less than
    /// ncols() rows).
    size_t nrows_cache_block () const { return nrows_cache_block_; }
  };

} // namespace TSQR


#endif // __TSQR_CacheBlocker_hpp
