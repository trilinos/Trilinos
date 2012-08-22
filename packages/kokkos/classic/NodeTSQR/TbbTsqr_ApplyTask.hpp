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

#ifndef __TSQR_TBB_ApplyTask_hpp
#define __TSQR_TBB_ApplyTask_hpp

#include <tbb/task.h>
#include <TbbTsqr_Partitioner.hpp>
#include <Tsqr_SequentialTsqr.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {
    
    /// \class ApplyTask
    /// \brief TBB task for recursive TSQR "apply Q factor" phase.
    ///
    template< class LocalOrdinal, class Scalar, class TimerType >
    class ApplyTask : public tbb::task {
    public:
      typedef MatView<LocalOrdinal, Scalar> mat_view;
      typedef ConstMatView<LocalOrdinal, Scalar> const_mat_view;
      typedef std::pair<mat_view, mat_view> split_t;
      typedef std::pair<const_mat_view, const_mat_view> const_split_t;
      typedef std::pair<const_mat_view, mat_view> top_blocks_t;
      typedef std::vector<top_blocks_t> array_top_blocks_t;

      /// \typedef SeqOutput
      /// Result of SequentialTsqr for each thread.
      typedef typename SequentialTsqr<LocalOrdinal, Scalar>::FactorOutput SeqOutput;
      /// \typedef ParOutput
      ///
      /// Array of ncores "local tau arrays" from parallel TSQR.
      /// (Local Q factors are stored in place.)
      typedef std::vector<std::vector<Scalar> > ParOutput;
      /// \typedef FactorOutput
      /// Result of SequentialTsqr for the data on each thread,
      /// and the result of combining the threads' data.
      typedef typename std::pair<std::vector<SeqOutput>, ParOutput> FactorOutput;

      /// \brief Constructor.
      ///
      /// \note The timing references are only modified by one thread
      ///   at a time; recursive calls use distinct references and
      ///   combine the results.
      ApplyTask (const size_t P_first__, 
		 const size_t P_last__,
		 ConstMatView<LocalOrdinal, Scalar> Q,
		 MatView<LocalOrdinal, Scalar> C,
		 array_top_blocks_t& top_blocks, 
		 const FactorOutput& factor_output,
		 const SequentialTsqr<LocalOrdinal, Scalar>& seq,
		 double& my_seq_timing,
		 double& min_seq_timing,
		 double& max_seq_timing,
		 const bool contiguous_cache_blocks) :
	P_first_ (P_first__), 
	P_last_ (P_last__), 
	Q_ (Q), 
	C_ (C),
	top_blocks_ (top_blocks), 
	factor_output_ (factor_output), 
	seq_ (seq), 
	apply_type_ (ApplyType::NoTranspose), // FIXME: modify to support Q^T and Q^H
	my_seq_timing_ (my_seq_timing),
	min_seq_timing_ (min_seq_timing),
	max_seq_timing_ (max_seq_timing),
	contiguous_cache_blocks_ (contiguous_cache_blocks)
      {}

      tbb::task* execute () 
      {
	if (P_first_ > P_last_ || Q_.empty() || C_.empty())
	  return NULL;
	else if (P_first_ == P_last_)
	  {
	    execute_base_case ();
	    return NULL;
	  }
	else
	  {
	    // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
	    const size_t P_mid = (P_first_ + P_last_) / 2;
	    const_split_t Q_split = 
	      partitioner_.split (Q_, P_first_, P_mid, P_last_,
				  contiguous_cache_blocks_);
	    split_t C_split = 
	      partitioner_.split (C_, P_first_, P_mid, P_last_,
				  contiguous_cache_blocks_);

	    // The partitioner may decide that the current blocks Q_
	    // and C_ have too few rows to be worth splitting.  In
	    // that case, Q_split.second and C_split.second (the
	    // bottom block) will be empty.  We can deal with this by
	    // treating it as the base case.
	    if (Q_split.second.empty() || Q_split.second.nrows() == 0)
	      {
		execute_base_case ();
		return NULL;
	      }

	    double top_timing;
	    double top_min_timing = 0.0;
	    double top_max_timing = 0.0;
	    double bot_timing;
	    double bot_min_timing = 0.0;
	    double bot_max_timing = 0.0;

	    apply_pair (P_first_, P_mid+1);
	    ApplyTask& topTask = *new( allocate_child() )
	      ApplyTask (P_first_, P_mid, Q_split.first, C_split.first,
			 top_blocks_, factor_output_, seq_, 
			 top_timing, top_min_timing, top_max_timing,
			 contiguous_cache_blocks_);
	    ApplyTask& botTask = *new( allocate_child() )
	      ApplyTask (P_mid+1, P_last_, Q_split.second, C_split.second,
			 top_blocks_, factor_output_, seq_,
			 bot_timing, bot_min_timing, bot_max_timing,
			 contiguous_cache_blocks_);

	    set_ref_count (3); // 3 children (2 + 1 for the wait)
	    spawn (topTask);
	    spawn_and_wait_for_all (botTask);

	    top_min_timing = (top_min_timing == 0.0) ? top_timing : top_min_timing;
	    top_max_timing = (top_max_timing == 0.0) ? top_timing : top_max_timing;

	    bot_min_timing = (bot_min_timing == 0.0) ? bot_timing : bot_min_timing;
	    bot_max_timing = (bot_max_timing == 0.0) ? bot_timing : bot_max_timing;

	    min_seq_timing_ = std::min (top_min_timing, bot_min_timing);
	    max_seq_timing_ = std::min (top_max_timing, bot_max_timing);

	    return NULL;
	  }
      }

    private:
      size_t P_first_, P_last_;
      const_mat_view Q_;
      mat_view C_;
      array_top_blocks_t& top_blocks_;
      const FactorOutput& factor_output_;
      SequentialTsqr<LocalOrdinal, Scalar> seq_;
      TSQR::ApplyType apply_type_;
      TSQR::Combine<LocalOrdinal, Scalar> combine_;
      Partitioner<LocalOrdinal, Scalar> partitioner_;
      double& my_seq_timing_;
      double& min_seq_timing_;
      double& max_seq_timing_;
      bool contiguous_cache_blocks_;

      void 
      execute_base_case ()
      {
	TimerType timer("");
	timer.start();
	const std::vector<SeqOutput>& seq_outputs = factor_output_.first;
	seq_.apply (apply_type_, Q_.nrows(), Q_.ncols(), 
		    Q_.get(), Q_.lda(), seq_outputs[P_first_], 
		    C_.ncols(), C_.get(), C_.lda(), 
		    contiguous_cache_blocks_);
	my_seq_timing_ = timer.stop();
      }

      void 
      apply_pair (const size_t P_top, 
		  const size_t P_bot) 
      {
	if (P_top == P_bot) 
	  throw std::logic_error("apply_pair: should never get here!");

	const_mat_view& Q_bot = top_blocks_[P_bot].first;
	mat_view& C_top = top_blocks_[P_top].second;
	mat_view& C_bot = top_blocks_[P_bot].second;

	const ParOutput& par_output = factor_output_.second;
	const std::vector<Scalar>& tau = par_output[P_bot];
	std::vector<Scalar> work (C_top.ncols());
	combine_.apply_pair (apply_type_, C_top.ncols(), Q_bot.ncols(), 
			     Q_bot.get(), Q_bot.lda(), &tau[0],
			     C_top.get(), C_top.lda(), 
			     C_bot.get(), C_bot.lda(), &work[0]);
      }

    };

  } // namespace TBB
} // namespace TSQR


#endif // __TSQR_TBB_ApplyTask_hpp
