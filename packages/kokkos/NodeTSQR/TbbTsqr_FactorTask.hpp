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

#ifndef __TSQR_TBB_FactorTask_hpp
#define __TSQR_TBB_FactorTask_hpp

#include <tbb/task.h>
#include <TbbTsqr_Partitioner.hpp>
#include <Tsqr_SequentialTsqr.hpp>
#include <Teuchos_TestForException.hpp>
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {
    
    /// \class FactorTask
    /// \brief TBB task for recursive TSQR factorization phase.
    ///
    template<class LocalOrdinal, class Scalar, class TimerType>
    class FactorTask : public tbb::task {
    public:
      typedef MatView<LocalOrdinal, Scalar> mat_view;
      typedef ConstMatView<LocalOrdinal, Scalar> const_mat_view;
      typedef std::pair<mat_view, mat_view> split_t;
      typedef std::pair<const_mat_view, const_mat_view> const_split_t;

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
      FactorTask (const size_t P_first__, 
		  const size_t P_last__,
		  mat_view A,
		  mat_view* const A_top_ptr,
		  std::vector<SeqOutput>& seq_outputs,
		  ParOutput& par_output,
		  const SequentialTsqr<LocalOrdinal, Scalar>& seq,
		  double& my_seq_timing,
		  double& min_seq_timing,
		  double& max_seq_timing,
		  const bool contiguous_cache_blocks) :
	P_first_ (P_first__),
	P_last_ (P_last__),
	A_ (A),
	A_top_ptr_ (A_top_ptr),
	seq_outputs_ (seq_outputs),
	par_output_ (par_output),
	seq_ (seq),
	contiguous_cache_blocks_ (contiguous_cache_blocks),
	my_seq_timing_ (my_seq_timing),
	min_seq_timing_ (min_seq_timing),
	max_seq_timing_ (max_seq_timing)
      {}

      tbb::task* execute () 
      {
	if (P_first_ > P_last_ || A_.empty())
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
	    split_t A_split = 
	      partitioner_.split (A_, P_first_, P_mid, P_last_,
				  contiguous_cache_blocks_);
	    // The partitioner may decide that the current block A_
	    // has too few rows to be worth splitting.  In that case,
	    // A_split.second (the bottom block) will be empty.  We
	    // can deal with this by treating it as the base case.
	    if (A_split.second.empty() || A_split.second.nrows() == 0)
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

	    FactorTask& topTask = *new( allocate_child() )
	      FactorTask (P_first_, P_mid, A_split.first, A_top_ptr_, 
			  seq_outputs_, par_output_, seq_,
			  top_timing, top_min_timing, top_max_timing,
			  contiguous_cache_blocks_);
	    // After the task finishes, A_bot will be set to the topmost
	    // partition of A_split.second.  This will let us combine
	    // the two subproblems (using factor_pair()) after their
	    // tasks complete.
	    mat_view A_bot;
	    FactorTask& botTask = *new( allocate_child() )
	      FactorTask (P_mid+1, P_last_, A_split.second, &A_bot, 
			  seq_outputs_, par_output_, seq_,
			  bot_timing, bot_min_timing, bot_max_timing,
			  contiguous_cache_blocks_);
	    set_ref_count (3); // 3 children (2 + 1 for the wait)
	    spawn (topTask);
	    spawn_and_wait_for_all (botTask);
	    
	    // Combine the two results
	    factor_pair (P_first_, P_mid+1, *A_top_ptr_, A_bot);

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
      const size_t P_first_, P_last_;
      mat_view A_;
      mat_view* const A_top_ptr_;
      std::vector<SeqOutput>& seq_outputs_;
      ParOutput& par_output_;
      SequentialTsqr<LocalOrdinal, Scalar> seq_;
      TSQR::Combine<LocalOrdinal, Scalar> combine_;
      Partitioner<LocalOrdinal, Scalar> partitioner_;
      const bool contiguous_cache_blocks_;
      double& my_seq_timing_;
      double& min_seq_timing_;
      double& max_seq_timing_;

      void 
      factor_pair (const size_t P_top,
		   const size_t P_bot,
		   mat_view& A_top, // different than A_top_
		   mat_view& A_bot)
      {
	const char thePrefix[] = "TSQR::TBB::Factor::factor_pair: ";
	TEST_FOR_EXCEPTION(P_top == P_bot, std::logic_error,
			   thePrefix << "Should never get here! P_top == P_bot (= " 
			   << P_top << "), that is, the indices of the thread "
			   "partitions are the same.");
	// We only read and write the upper ncols x ncols triangle of
	// each block.
	TEST_FOR_EXCEPTION(A_top.ncols() != A_bot.ncols(), std::logic_error,
			   thePrefix << "The top cache block A_top is " 
			   << A_top.nrows() << " x " << A_top.ncols() 
			   << ", and the bottom cache block A_bot is "
			   << A_bot.nrows() << " x " << A_bot.ncols() 
			   << "; this means we can't factor [A_top; A_bot].");
	const LocalOrdinal ncols = A_top.ncols();
	std::vector<Scalar>& tau = par_output_[P_bot];
	std::vector<Scalar> work (ncols);
	combine_.factor_pair (ncols, A_top.get(), A_top.lda(),
			      A_bot.get(), A_bot.lda(), &tau[0], &work[0]);
      }

      void
      execute_base_case () 
      {
	TimerType timer("");
	timer.start();
	seq_outputs_[P_first_] = 
	  seq_.factor (A_.nrows(), A_.ncols(), A_.get(), 
		       A_.lda(), contiguous_cache_blocks_);
	// Assign the topmost cache block of the current partition to
	// *A_top_ptr_.  Every base case invocation does this, so that
	// we can combine subproblems.  The root task also does this,
	// but for a different reason: so that we can extract the R
	// factor, once we're done with the factorization.
	*A_top_ptr_ = seq_.top_block (A_, contiguous_cache_blocks_);
	my_seq_timing_ = timer.stop();
      }
    };
  } // namespace TBB
} // namespace TSQR

#endif // __TSQR_TBB_FactorTask_hpp
