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

#ifndef __TSQR_TBB_FillWithZerosTask_hpp
#define __TSQR_TBB_FillWithZerosTask_hpp

#include <tbb/task.h>
#include <TbbTsqr_Partitioner.hpp>
#include <Tsqr_SequentialTsqr.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {

    /// \class FillWithZerosTask
    /// \brief TBB task for recursive TSQR "fill with zeros" phase.
    ///
    template<class LocalOrdinal, class Scalar>
    class FillWithZerosTask : public tbb::task {
    private:
      typedef MatView<LocalOrdinal, Scalar> mat_view;
      typedef std::pair<mat_view, mat_view> split_type;

    public:
      FillWithZerosTask (const size_t P_first, 
			 const size_t P_last,
			 MatView<LocalOrdinal, Scalar> C,
			 const SequentialTsqr<LocalOrdinal, Scalar>& seq,
			 const bool contiguous_cache_blocks = false)
	: P_first_ (P_first), 
	  P_last_ (P_last), 
	  C_ (C),
	  seq_ (seq), 
	  contiguous_cache_blocks_ (contiguous_cache_blocks)
      {}

      tbb::task* execute () 
      {
	if (P_first_ > P_last_ || C_.empty())
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
	    split_type C_split = 
	      partitioner_.split (C_, P_first_, P_mid, P_last_, 
				  contiguous_cache_blocks_);
	    // The partitioner may decide that the current block C_
	    // has too few rows to be worth splitting.  In that case,
	    // C_split.second (the bottom block) will be empty.  We
	    // can deal with this by treating it as the base case.
	    if (C_split.second.empty() || C_split.second.nrows() == 0)
	      {
		execute_base_case ();
		return NULL;
	      }

	    // "c": continuation task
	    tbb::empty_task& c = 
	      *new( allocate_continuation() ) tbb::empty_task;
	    // Recurse on the split
	    FillWithZerosTask& topTask = *new( c.allocate_child() )
	      FillWithZerosTask (P_first_, P_mid, C_split.first, seq_, 
				 contiguous_cache_blocks_);
	    FillWithZerosTask& botTask = *new( c.allocate_child() )
	      FillWithZerosTask (P_mid+1, P_last_, C_split.second, seq_, 
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
      mat_view C_;
      SequentialTsqr<LocalOrdinal, Scalar> seq_;
      Partitioner<LocalOrdinal, Scalar> partitioner_;
      bool contiguous_cache_blocks_;

      void 
      execute_base_case ()
      {
	// Fill my partition with zeros.
	seq_.fill_with_zeros (C_.nrows(), C_.ncols(), C_.get(), 
			      C_.lda(), contiguous_cache_blocks_);
      }
    };
  } // namespace TBB
} // namespace TSQR


#endif // __TSQR_TBB_FillWithZerosTask_hpp
