// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef __TSQR_TBB_RevealRankTask_hpp
#define __TSQR_TBB_RevealRankTask_hpp

#include <tbb/task.h>
#include <TbbTsqr_Partitioner.hpp>
#include <Tsqr_SequentialTsqr.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {
    
    template< class LocalOrdinal, class Scalar >
    class RevealRankTask : public tbb::task {
    public:
      typedef MatView< LocalOrdinal, Scalar > mat_view;
      typedef ConstMatView< LocalOrdinal, Scalar > const_mat_view;
      typedef std::pair< mat_view, mat_view > split_type;
      typedef SequentialTsqr< LocalOrdinal, Scalar > seq_tsqr_type;

      RevealRankTask (const size_t P_first, 
		      const size_t P_last,
		      const mat_view& Q,
		      const const_mat_view& U,
		      const seq_tsqr_type& seq,
		      const bool contiguous_cache_blocks) :
	P_first_ (P_first), 
	P_last_ (P_last), 
	Q_ (Q),
	U_ (U),
	seq_ (seq),
	contiguous_cache_blocks_ (contiguous_cache_blocks)
      {}

      tbb::task* execute () {
	using tbb::task;

	if (P_first_ > P_last_ || Q_.empty())
	  return NULL; // shouldn't get here, but just in case...
	else if (P_first_ == P_last_)
	  {
	    // Use SequentialTsqr to compute Q*U for this core's local
	    // part of Q.  The method is called "Q_times_B" so that it
	    // doesn't suggest any orthogonality of the B input
	    // matrix, though in this case B is U and U is orthogonal.
	    seq_.Q_times_B (Q_.nrows(), Q_.ncols(), Q_.get(), Q_.lda(),
			    U_.get(), U_.lda(), contiguous_cache_blocks_);
	    return NULL;
	  }
	else
	  {
	    // "c": continuation task
	    tbb::empty_task& c = *new( allocate_continuation() ) tbb::empty_task;

	    // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
	    const size_t P_mid = (P_first_ + P_last_) / 2;
	    split_type out_split = 
	      partitioner_.split (Q_, P_first_, P_mid, P_last_, 
				  contiguous_cache_blocks_);

	    RevealRankTask& topTask = *new( c.allocate_child() )
	      RevealRankTask (P_first_, P_mid, out_split.first, seq_, 
			      contiguous_cache_blocks_);
	    RevealRankTask& botTask = *new( c.allocate_child() )
	      RevealRankTask (P_mid+1, P_last_, out_split.second, seq_,
			      contiguous_cache_blocks_);
	    // Set reference count of parent (in this case, the
	    // continuation task) to 2 (since 2 children -- no
	    // additional task since no waiting).
	    c.set_ref_count (2);
	    c.spawn (botTask);
	    return &topTask; // scheduler bypass optimization
	  }
      }

    private:
      size_t P_first_, P_last_;
      mat_view Q_;
      const_mat_view U_;
      SequentialTsqr< LocalOrdinal, Scalar > seq_;
      Partitioner< LocalOrdinal, Scalar > partitioner_;
      bool contiguous_cache_blocks_;
    };

  } // namespace TBB
} // namespace TSQR


#endif // __TSQR_TBB_RevealRankTask_hpp
