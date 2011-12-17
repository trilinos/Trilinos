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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_TBB_Partitioner_hpp
#define __TSQR_TBB_Partitioner_hpp

#include <Tsqr_MatView.hpp>

#include <cstring> // size_t
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {

    template< class Ordinal, class Scalar >
    class Partitioner {
    private:
      bool
      should_split (const Ordinal nrows,
		    const Ordinal ncols,
		    const size_t num_partitions) const
      {
	using std::invalid_argument;
	using std::ostringstream;

	if (nrows < ncols)
	  {
	    ostringstream os;
	    os << "Partitioner::should_split: nrows (= " << nrows 
	       << ") < ncols (= " << ncols << ")";
	    throw invalid_argument (os.str());
	  }
	else if (num_partitions == 0)
	  {
	    ostringstream os;
	    os << "Partitioner::should_split: nrows (= " << nrows 
	       << ") < ncols (= " << ncols << ")";
	    throw invalid_argument (os.str());
	  }
	// FIXME (mfh 11 Jul 2010) Need more overflow checks here.
	return static_cast<size_t>(nrows) / num_partitions >= static_cast<size_t>(ncols);
      }	

    public:
      /// Partition into [P_first, P_mid] and [P_mid+1, P_last].  The
      /// base case is reached when the second returned MatrixViewType
      /// is empty.
      template< class MatrixViewType >
      std::pair< MatrixViewType, MatrixViewType >
      split (const MatrixViewType& A,
	     const size_t P_first,
	     const size_t P_mid,
	     const size_t P_last,
	     const bool contiguous_cache_blocks) const
      {
	typedef typename MatrixViewType::ordinal_type ordinal_type;
	typedef typename MatrixViewType::pointer_type pointer_type;

	const size_t num_partitions_top = P_mid - P_first + 1;
	//const size_t num_partitions_bottom = P_last - P_mid;
	const size_t num_partitions = P_last - P_first + 1;
	const ordinal_type nrows = A.nrows();
	const ordinal_type ncols = A.ncols();
	
	if (! should_split (nrows, ncols, num_partitions))
	  return std::make_pair (MatrixViewType(A), MatrixViewType());
	else
	  {
	    const ordinal_type num_rows_partition = nrows / num_partitions;
	    const ordinal_type remainder = nrows % num_partitions;
	    
	    // Top partition gets the remainder rows.  Doing the
	    // multiplication before the division might make it more
	    // likely to avoid truncating the fraction, but may cause
	    // overflow of ordinal_type.  
	    const ordinal_type num_rows_top = 
	      num_rows_partition * num_partitions_top + remainder;
	    const ordinal_type num_rows_bot = nrows - num_rows_top;

	    // We don't call (Const)MatView::split_top(), because that
	    // is for splitting off a single cache block.  Each half
	    // of the split may contain more than one cache block.
	    if (contiguous_cache_blocks)
	      {
		pointer_type A_bot_ptr = A.get() + num_rows_top * ncols;
		MatrixViewType A_top (num_rows_top, ncols, A.get(), num_rows_top);
		MatrixViewType A_bot (num_rows_bot, ncols, A_bot_ptr, num_rows_bot);
		return std::make_pair (A_top, A_bot);
	      }
	    else
	      {
		pointer_type A_bot_ptr = A.get() + num_rows_top;
		MatrixViewType A_top (num_rows_top, ncols, A.get(), A.lda());
		MatrixViewType A_bot (num_rows_bot, ncols, A_bot_ptr, A.lda());
		return std::make_pair (A_top, A_bot);
	      }
	  }
      }
    }; // class Partitioner
  } // namespace TBB
} // namespace TSQR

#endif // __TSQR_TBB_Partitioner_hpp
