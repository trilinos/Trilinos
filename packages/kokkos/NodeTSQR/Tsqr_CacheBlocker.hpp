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

// #include <iostream>
#include <sstream>
#include <stdexcept>

// #define CACHE_BLOCKER_DEBUG 1
// #ifdef CACHE_BLOCKER_DEBUG
// #  undef CACHE_BLOCKER_DEBUG
// #endif // CACHE_BLOCKER_DEBUG
#ifdef CACHE_BLOCKER_DEBUG
#  include <iostream>
using std::cerr;
using std::endl;
#endif // CACHE_BLOCKER_DEBUG

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class CacheBlocker
  /// 
  /// A CacheBlocker uses a particular cache blocking strategy ("class
  /// CBS") to know how to separate an nrows by ncols matrix into
  /// cache blocks, which in this case are row blocks.  We can use
  /// this strategy to cache block any matrix of nrows rows, using the
  /// same distribution of rows to blocks, regardless of whether that
  /// latter matrix has a different number of columns than ncols.
  template< class Ordinal, class Scalar > // class CBS=CacheBlockingStrategy< Ordinal, Scalar > 
  class CacheBlocker {
  private:
    typedef MatView< Ordinal, Scalar > mat_view;
    typedef ConstMatView< Ordinal, Scalar > const_mat_view;

    void
    validate () 
    {
      if (nrows_cache_block_ < ncols_)
	{
	  std::ostringstream os;
	  os << "Cache block is too small: only " << nrows_cache_block_
	     << " rows fit, but need at least as many rows as columns (= "
	     << ncols_ << ")";
	  throw std::logic_error (os.str());
	}
    }

  public:
    /// Constructor
    ///
    /// \param num_rows Number of rows in the matrix corresponding to
    ///   the cache blocking scheme
    /// \param num_cols Number of columns in the matrix corresponding
    ///   to the cache blocking scheme 
    /// \param strategy Cache blocking strategy object
    CacheBlocker (const Ordinal num_rows,
		  const Ordinal num_cols,
		  const CacheBlockingStrategy< Ordinal, Scalar >& strategy) :
      nrows_ (num_rows), 
      ncols_ (num_cols), 
      strategy_ (strategy)
    {
#ifdef CACHE_BLOCKER_DEBUG
      cerr << "CacheBlocker:" << endl
	   << "# rows = " << nrows() << endl
	   << "# cols = " << ncols() << endl
	   << "Cache block size (bytes) = " << strategy_.cache_block_size()
	   << endl << endl;
#endif
      const Ordinal temp_nrows_cache_block = 
	strategy_.cache_block_num_rows (ncols());
#ifdef CACHE_BLOCKER_DEBUG
      cerr << "Strategy says: # rows per cache block = "
	   << temp_nrows_cache_block << endl;
#endif
      nrows_cache_block_ = temp_nrows_cache_block;
      validate ();
    }

    /// Copy constructor
    CacheBlocker (const CacheBlocker& rhs) :
      nrows_ (rhs.nrows()), 
      ncols_ (rhs.ncols()), 
      strategy_ (rhs.strategy_)
    {}

    /// Assignment operator
    CacheBlocker& operator= (const CacheBlocker& rhs) {
      nrows_ = rhs.nrows();
      ncols_ = rhs.ncols();
      strategy_ = rhs.strategy_;
      return *this;
    }

    /// Cache block size in bytes
    size_t cache_block_size () const { return strategy_.cache_block_size(); }

    /// Number of rows in the matrix corresponding to the cache blocking scheme
    Ordinal nrows () const { return nrows_; }

    /// Number of columns in the matrix corresponding to the cache blocking scheme
    Ordinal ncols () const { return ncols_; }

    /// Return the topmost cache block of A, where the number of rows
    /// in each cache block is chosen according to the number of
    /// columns with which this CacheBlocker was set up (rather than
    /// the number of columns in A, which may not be the same).  The
    /// idea is to have the number and distribution of rows in the
    /// cache blocks be the same as the original nrows() by ncols()
    /// matrix with which this CacheBlocker was initialized.
    template< class MatrixViewType >
    MatrixViewType
    split_top_block (MatrixViewType& A, const bool contiguous_cache_blocks) const
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      const ordinal_type nrows_top = 
	strategy_.top_block_split_nrows (A.nrows(), ncols(), 
					 nrows_cache_block());
      // split_top() modifies A
      return A.split_top (nrows_top, contiguous_cache_blocks);
    }

    /// Return the topmost cache block of A.  A is copied so it won't
    /// be modified.
    template< class MatrixViewType >
    MatrixViewType
    top_block (const MatrixViewType& A, const bool contiguous_cache_blocks) const
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      const ordinal_type nrows_top = 
	strategy_.top_block_split_nrows (A.nrows(), ncols(), 
					 nrows_cache_block());
      MatrixViewType A_copy (A);
      return A_copy.split_top (nrows_top, contiguous_cache_blocks);
    }

    template< class MatrixViewType >
    MatrixViewType
    split_bottom_block (MatrixViewType& A, const bool contiguous_cache_blocks) const
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      const ordinal_type nrows_bottom = 
	strategy_.bottom_block_split_nrows (A.nrows(), ncols(), 
					    nrows_cache_block());
      // split_bottom() modifies A
      return A.split_bottom (nrows_bottom, contiguous_cache_blocks);
    }

    /// Fill the entries of A with Scalar(0).
    ///
    /// \note This method only works if MatrixViewType has a method
    ///   "fill(const Scalar)" or "fill(const Scalar&)".  The
    ///   intention is that the method be non-const and that it fill
    ///   in the entries of the matrix with Scalar(0), though
    ///   syntactically the method could even be const and / or do
    ///   something else entirely.
    template< class MatrixViewType >
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
	  // This call modifies A, but that's OK since we passed the
	  // input view by copy, not by reference.
	  MatrixViewType A_cur = split_top_block (A, contiguous_cache_blocks);
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

    /// "Un"-cache-block the given A_in matrix, writing the results to
    /// A_out.
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


    void
    fill_with_zeros (const Ordinal num_rows,
		     const Ordinal num_cols,
		     Scalar A[],
		     const Ordinal lda, 
		     const bool contiguous_cache_blocks)
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


  private:
    Ordinal nrows_;
    Ordinal ncols_;
    CacheBlockingStrategy< Ordinal, Scalar > strategy_;
    Ordinal nrows_cache_block_;

    /// Number of rows in a typical cache block (not every cache block
    /// has this many rows -- some may have more (up to but not
    /// including nrows_cache_block() + ncols()), and some may have
    /// less (whatever is left over at the bottom of the matrix, but
    /// no less than ncols() rows).
    size_t nrows_cache_block () const { return nrows_cache_block_; }
  };

} // namespace TSQR


#endif // __TSQR_CacheBlocker_hpp
