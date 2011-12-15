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

#ifndef __TSQR_TBB_UnCacheBlockTask_hpp
#define __TSQR_TBB_UnCacheBlockTask_hpp

#include <tbb/task.h>
#include <TbbTsqr_Partitioner.hpp>
#include <Tsqr_SequentialTsqr.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {
    
    /// \class UnCacheBlockTask
    /// \brief TBB task for recursive TSQR un-(cache blocking) phase.
    /// 
    /// "Un-(cache blocking)" here means copying the input matrix,
    /// which is stored with contiguous cache blocks, to the output
    /// matrix, which is stored with noncontiguous cache blocks.
    ///
    template<class LocalOrdinal, class Scalar>
    class UnCacheBlockTask : public tbb::task {
    public:
      typedef MatView< LocalOrdinal, Scalar > mat_view;
      typedef ConstMatView< LocalOrdinal, Scalar > const_mat_view;
      typedef std::pair< mat_view, mat_view > split_t;
      typedef std::pair< const_mat_view, const_mat_view > const_split_t;

      UnCacheBlockTask (const size_t P_first__, 
			const size_t P_last__,
			mat_view& A_out,
			const_mat_view& A_in,
			const SequentialTsqr<LocalOrdinal, Scalar>& seq) :
	P_first_ (P_first__), 
	P_last_ (P_last__), 
	A_out_ (A_out), 
	A_in_ (A_in), 
	seq_ (seq)
      {}

      tbb::task* execute () 
      {
	using tbb::task;

	if (P_first_ > P_last_ || A_out_.empty() || A_in_.empty())
	  return NULL;
	else if (P_first_ == P_last_)
	  {
	    execute_base_case ();
	    return NULL;
	  }
	else
	  {
	    // Recurse on two intervals: [P_first, P_mid] and
	    // [P_mid+1, P_last].
	    const size_t P_mid = (P_first_ + P_last_) / 2;
	    split_t out_split = 
	      partitioner_.split (A_out_, P_first_, P_mid, P_last_, false);
	    const_split_t in_split = 
	      partitioner_.split (A_in_, P_first_, P_mid, P_last_, true);

	    // The partitioner may decide that the current blocks
	    // A_out_ and A_in_ have too few rows to be worth
	    // splitting.  (It should split both A_out_ and A_in_ in
	    // the same way.)  In that case, out_split.second and
	    // in_split.second (the bottom block) will be empty.  We
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
	    UnCacheBlockTask& topTask = *new( c.allocate_child() )
	      UnCacheBlockTask (P_first_, P_mid, out_split.first, 
			      in_split.first, seq_);
	    UnCacheBlockTask& botTask = *new( c.allocate_child() )
	      UnCacheBlockTask (P_mid+1, P_last_, out_split.second, 
			      in_split.second, seq_);
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
      mat_view A_out_;
      const_mat_view A_in_;
      SequentialTsqr<LocalOrdinal, Scalar> seq_;
      Partitioner<LocalOrdinal, Scalar> partitioner_;

      void
      execute_base_case ()
      {
	seq_.un_cache_block (A_out_.nrows(), A_out_.ncols(), 
			     A_out_.get(), A_out_.lda(), A_in_.get());
      }
    };

  } // namespace TBB
} // namespace TSQR


#endif // __TSQR_TBB_UnCacheBlockTask_hpp
