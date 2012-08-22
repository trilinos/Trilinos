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

#ifndef __TSQR_TBB_RevealRankTask_hpp
#define __TSQR_TBB_RevealRankTask_hpp

#include <tbb/task.h>
#include <TbbTsqr_Partitioner.hpp>
#include <Tsqr_SequentialTsqr.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {

    /// \class RevealRankTask
    /// \brief TBB task for recursive TSQR "rank-revealing" phase.
    ///
    /// This part of the factorization doesn't actually reveal the
    /// rank in parallel; we assume that this has already been done
    /// and the columns of U form a basis for the column space of the
    /// R factor (in the QR factorization of the original matrix).
    /// All we need to do here is compute Q*U in parallel, respecting
    /// the original partitioning and cache blocking scheme.
    template<class LocalOrdinal, class Scalar>
    class RevealRankTask : public tbb::task {
    public:
      typedef MatView<LocalOrdinal, Scalar> mat_view;
      typedef ConstMatView<LocalOrdinal, Scalar> const_mat_view;
      typedef std::pair<mat_view, mat_view> split_type;
      typedef SequentialTsqr<LocalOrdinal, Scalar> seq_tsqr_type;

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

      void 
      execute_base_case ()
      {
	// Use SequentialTsqr to compute Q*U for this core's local
	// part of Q.  The method is called "Q_times_B" so that it
	// doesn't suggest any orthogonality of the B input matrix,
	// though in this case B is U and U is orthogonal
	// (resp. unitary if Scalar is complex).
	seq_.Q_times_B (Q_.nrows(), Q_.ncols(), Q_.get(), Q_.lda(),
			U_.get(), U_.lda(), contiguous_cache_blocks_);
      }

      tbb::task* execute () 
      {
	using tbb::task;

	if (P_first_ > P_last_ || Q_.empty())
	  return NULL; // shouldn't get here, but just in case...
	else if (P_first_ == P_last_)
	  {
	    execute_base_case ();
	    return NULL;
	  }
	else
	  {
	    // Recurse on two intervals: [P_first, P_mid] and
	    // [P_mid+1, P_last]
	    const size_t P_mid = (P_first_ + P_last_) / 2;
	    split_type out_split = 
	      partitioner_.split (Q_, P_first_, P_mid, P_last_, 
				  contiguous_cache_blocks_);
	    // The partitioner may decide that the current block Q_
	    // has too few rows to be worth splitting.  In that case,
	    // out_split.second (the bottom block) will be empty.  We
	    // can deal with this by treating it as the base case.
	    if (out_split.second.empty() || out_split.second.nrows() == 0)
	      {
		execute_base_case ();
		return NULL;
	      }

	    // "c": continuation task
	    tbb::empty_task& c = 
	      *new( allocate_continuation() ) tbb::empty_task;
	    // Recurse on the split
	    RevealRankTask& topTask = *new( c.allocate_child() )
	      RevealRankTask (P_first_, P_mid, out_split.first, U_, 
			      seq_, contiguous_cache_blocks_);
	    RevealRankTask& botTask = *new( c.allocate_child() )
	      RevealRankTask (P_mid+1, P_last_, out_split.second, U_,
			      seq_, contiguous_cache_blocks_);
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
      SequentialTsqr<LocalOrdinal, Scalar> seq_;
      Partitioner<LocalOrdinal, Scalar> partitioner_;
      bool contiguous_cache_blocks_;
    };

  } // namespace TBB
} // namespace TSQR


#endif // __TSQR_TBB_RevealRankTask_hpp
