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

#ifndef __TSQR_CacheBlockingStrategy_hpp
#define __TSQR_CacheBlockingStrategy_hpp

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility> // std::pair

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class CacheBlockingStrategy
  ///
  /// Strategy object that helps CacheBlocker decide how to block up a
  /// given multivector.
  template< class LocalOrdinal, class Scalar >
  class CacheBlockingStrategy {
  public:
    /// Constructor
    ///
    /// \param cacheBlockSize [in] Cache block size in bytes.  If
    ///   zero, the cache block size is set to a reasonable default.
    CacheBlockingStrategy (const size_t cacheBlockSize = 0) :
      cache_block_size_ (default_cache_block_size (cacheBlockSize))
    {}

    ///
    /// Copy constructor
    ///
    CacheBlockingStrategy (const CacheBlockingStrategy& rhs) :
      cache_block_size_ (rhs.cache_block_size()) 
    {}

    ///
    /// Assignment operator
    ///
    CacheBlockingStrategy& operator= (const CacheBlockingStrategy& rhs) {
      cache_block_size_ = rhs.cache_block_size();
      return *this;
    }

    ///
    /// Return the cache block size in bytes
    ///
    size_t cache_block_size () const { return cache_block_size_; }

    /// Return true if the two strategies are the same
    ///
    /// \note Currently this means only that they use the same cache
    ///   block size.
    bool operator== (const CacheBlockingStrategy& rhs) const {
      return (cache_block_size() == rhs.cache_block_size());
    }

    /// Return true if the two strategies are not the same
    ///
    /// \note Currently this means only that they use different cache
    ///   block sizes.
    bool operator!= (const CacheBlockingStrategy& rhs) const {
      return (cache_block_size() != rhs.cache_block_size());
    }

    /// \return Row start index (zero-based) of my cache block, and
    ///   the number of rows in my cache block.  If the index is out
    ///   of range, the row start index will be nrows and the number
    ///   of rows will be zero.
    std::pair< LocalOrdinal, LocalOrdinal >
    cache_block (const LocalOrdinal index,
		 const LocalOrdinal nrows,
		 const LocalOrdinal ncols,
		 const LocalOrdinal nrows_cache_block) const
    {
      const LocalOrdinal nblocks = nrows / nrows_cache_block;
      const LocalOrdinal remainder = nrows - nblocks * nrows_cache_block;
      LocalOrdinal my_row_start, my_nrows;

      my_row_start = index * nrows_cache_block;
      if (index < nblocks - 1)
	my_nrows = nrows_cache_block;
      else if (index == nblocks - 1)
	{
	  if (remainder > 0 && remainder < ncols)
	    my_nrows = nrows_cache_block + remainder;
	  else
	    my_nrows = nrows_cache_block;
	}
      else if (index == nblocks)
	{
	  if (remainder >= ncols)
	    my_nrows = remainder;
	  else 
	    my_nrows = 0;
	}
      else 
	my_nrows = 0;
      
      return std::make_pair (my_row_start, my_nrows);
    }

    /// Number of cache blocks into which to break up an nrows by
    /// ncols matrix, where the suggested number of rows per cache
    /// block is nrows_cache_block, and where no cache block has fewer
    /// than ncols rows.
    LocalOrdinal 
    num_cache_blocks (const LocalOrdinal nrows,
		      const LocalOrdinal ncols,
		      const LocalOrdinal nrows_cache_block) const
    {
      const LocalOrdinal quotient = nrows / nrows_cache_block; 
      const LocalOrdinal remainder = nrows - nrows_cache_block * quotient;
      const LocalOrdinal nblocks = (0 < remainder && remainder < ncols) ? (quotient+1) : quotient;

      return nblocks;
    }

    /// If we partition the nrows by ncols matrix A into [A_top;
    /// A_bot] with A_top being a cache block and A_bot being the rest
    /// of the matrix, return the number of rows that A_top should
    /// have.
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


    /// If we partition the nrows by ncols matrix A into [A_top;
    /// A_bot] with A_bot being a cache block and A_top being the rest
    /// of the matrix, return the number of rows that A_bot should
    /// have.
    /// 
    /// \return # of rows in top cache block A_bot
    LocalOrdinal
    bottom_block_split_nrows (const LocalOrdinal nrows,
			      const LocalOrdinal ncols,
			      const LocalOrdinal nrows_cache_block) const
    {
      // We want to split off the bottom block using the same
      // splitting as if we had split off as many top blocks of
      // nrows_cache_block rows as permissible.  The last block may
      // have fewer than nrows_cache_block rows, but it may not have
      // fewer than ncols rows (since we don't want any cache block to
      // have fewer rows than columns).

      if (nrows < nrows_cache_block + ncols)
	{
	  // The nrows by ncols matrix A is just big enough for one
	  // cache block.  Partition A == [A_top; A_bot] with A_top
	  // empty and A_bot == A.
	  return nrows;
	}
      else
	{
	  const LocalOrdinal quotient = nrows / nrows_cache_block;
	  const LocalOrdinal remainder = nrows - quotient * nrows_cache_block;

	  LocalOrdinal nrows_bottom;
	  if (remainder == 0 || remainder < ncols)
	    nrows_bottom = nrows_cache_block + remainder;
	  else if (remainder >= ncols)
	    nrows_bottom = remainder;
	  else
	    throw std::logic_error("Should never get here!");
	  return nrows_bottom;
	}
    }


    /// Return the cache block size in bytes, based on the given
    /// suggested size in bytes, and the size of the Scalar type.  If
    /// the input is zero, return a reasonable default size.
    ///
    /// \return Either the given suggested cache block size (in bytes),
    /// or a default cache block size (in bytes) (if the given
    /// suggested size == 0).
    size_t 
    default_cache_block_size (const size_t suggested_cache_size) const
    {
      if (suggested_cache_size == 0)
	{
	  // 64 KB: reasonable default upper bound for L2 cache
	  // capacity.  
	  const size_t default_cache_size = 65536;
	  const size_t default_nwords = default_cache_size / sizeof(Scalar);

	  // We want to be able to work on two blocks of size at least
	  // 4x4 in cache.  Otherwise TSQR is more or less pointless.
	  // If we can't do that, bail out; it's likely that Scalar is
	  // the wrong data type for TSQR.
	  if (default_nwords < 32)
	    {
	      std::ostringstream os;
	      os << "Error: sizeof(Scalar) == " << sizeof(Scalar) 
		 << ", which is too large for us to pick a reasonable "
		 << "default cache size.";
	      throw std::range_error (os.str());
	    }
	  else
	    return default_cache_size;
	}
      else
	{
	  const size_t nwords = suggested_cache_size / sizeof(Scalar);

	  // We want to be able to work on two blocks of size at least
	  // 4x4 in cache.  Otherwise TSQR is more or less pointless.
	  // If we can't do that, bail out; it's likely that Scalar is
	  // the wrong data type for TSQR.
	  if (nwords < 32)
	    {
	      std::ostringstream os;
	      os << "Error: suggested cache size " << suggested_cache_size
		 << " bytes is too small for sizeof(Scalar) == " 
		 << sizeof(Scalar) << " bytes";
	      throw std::invalid_argument (os.str());
	    }
	  else
	    return suggested_cache_size;
	}
    }


    /// Number of rows that a cache block should occupy, given a
    /// particular number of columns in the matrix to factor.  That
    /// cache block size is used to fix the number of rows per cache
    /// block also when applying Q to a matrix C.  We choose the cache
    /// block size so that when Q and C have the same number of
    /// columns, two cache blocks (one of Q and the other of C) will
    /// fit in cache.
    ///
    /// \param ncols [in] Number of columns in the matrix whose QR
    ///   factorization we compute using SequentialTsqr.
    LocalOrdinal 
    cache_block_num_rows (const LocalOrdinal ncols) const
    {
      // Suppose the cache can hold W words (of size sizeof(Scalar)
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
      // we use this formula to pick nrows_cache_block.

      const size_t cache_block_nwords = cache_block_size() / sizeof(Scalar);

      // Compute everything in size_t first, and cast to LocalOrdinal
      // at the end.  This may avoid overflow if the cache is very
      // large and/or LocalOrdinal is very small (say a short int).
      //
      // Also, since size_t is unsigned, make sure that the
      // subtractions don't make it negative.  If it does, then either
      // ncols is too big or the cache is too small.

      const size_t term1 = cache_block_nwords / (2*ncols);
      const size_t term2 = ncols / 2 + 1;
      if (term1 < term2)
	{
	  std::ostringstream os;
	  os << "Error:  While deciding on the number of rows in a cache "
	    "block for sequential TSQR, the specified number of columns " 
	     << ncols << " in the matrix to factor is too big for the "
	    "specified cache size, which can hold " << cache_block_nwords 
	     << " words of size sizeof(Scalar) == " << sizeof(Scalar) 
	     << " each.";
	  throw std::invalid_argument (os.str());
	}
      else 
	{
	  // The compiler can easily prove that term1 - term2 >= 0,
	  // since we've gotten to this point.  Of course that's
	  // assuming that C++ compilers are smart...
	  const size_t nrows_cache_block = term1 - term2;

	  // Make sure that nrows_cache_block fits in a LocalOrdinal type.
	  if (static_cast<size_t>(static_cast< LocalOrdinal >(nrows_cache_block)) != nrows_cache_block)
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
	    return static_cast< LocalOrdinal > (nrows_cache_block);
	}
    }

  private:
    /// Size in bytes of each cache block
    size_t cache_block_size_;
  };
} // namespace TSQR

#endif // __TSQR_CacheBlockingStrategy_hpp
