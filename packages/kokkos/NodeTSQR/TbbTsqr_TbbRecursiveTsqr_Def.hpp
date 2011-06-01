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

#ifndef __TSQR_TBB_TbbRecursiveTsqr_Def_hpp
#define __TSQR_TBB_TbbRecursiveTsqr_Def_hpp

#include <TbbTsqr_TbbRecursiveTsqr.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_Util.hpp>

// #define TBB_DEBUG 1
#ifdef TBB_DEBUG
#  include <iostream>
using std::cerr;
using std::endl;
#endif // TBB_DEBUG

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {

    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    explicit_Q_helper (const size_t P_first, 
		       const size_t P_last,
		       MatView< LocalOrdinal, Scalar >& Q_out,
		       const bool contiguous_cache_blocks) const
    {
      if (P_first > P_last || Q_out.empty())
	return;
      else if (P_first == P_last)
	{
	  CacheBlocker< LocalOrdinal, Scalar > 
	    blocker (Q_out.nrows(), Q_out.ncols(),
		     seq_.cache_blocking_strategy());
#ifdef TBB_DEBUG
	  cerr << "explicit_Q_helper: On P_first = " << P_first 
	       << ", filling Q_out with zeros:" << endl
	       << "Q_out is " << Q_out.nrows() << " x " << Q_out.ncols() 
	       << " with leading dimension " << Q_out.lda() << endl;
#endif // TBB_DEBUG
	  // Fill my partition with zeros.
	  blocker.fill_with_zeros (Q_out, contiguous_cache_blocks);

	  // If our partition is the first (topmost), fill it with
	  // the first Q_out.ncols() columns of the identity matrix.
	  if (P_first == 0)
	    {
	      // Fetch the topmost cache block of my partition.  Its
	      // leading dimension should be set correctly by
	      // top_block().
	      mat_view Q_out_top = 
		blocker.top_block (Q_out, contiguous_cache_blocks);

	      for (LocalOrdinal j = 0; j < Q_out_top.ncols(); ++j)
		Q_out_top(j,j) = Scalar(1);
	    }
	}
      else
	{
	  // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
	  const size_t P_mid = (P_first + P_last) / 2;
	  split_t Q_out_split =
	    partitioner_.split (Q_out, P_first, P_mid, P_last,
				contiguous_cache_blocks);
	  explicit_Q_helper (P_first, P_mid, Q_out_split.first, 
			     contiguous_cache_blocks);
	  explicit_Q_helper (P_mid+1, P_last, Q_out_split.second, 
			     contiguous_cache_blocks);
	}
    }


    template< class LocalOrdinal, class Scalar >
    MatView< LocalOrdinal, Scalar >
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    factor_helper (const size_t P_first, 
		   const size_t P_last,
		   const size_t depth,
		   MatView< LocalOrdinal, Scalar > A,
		   std::vector< typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::SeqOutput >& seq_outputs,
		   typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::ParOutput& par_outputs,
		   Scalar R[],
		   const LocalOrdinal ldr,
		   const bool contiguous_cache_blocks) const
    {
      mat_view A_top;
      if (P_first > P_last || A.empty())
	return A;
      else if (P_first == P_last)
	{
	  std::pair< SeqOutput, MatView< LocalOrdinal, Scalar > > results = 
	    seq_.factor (A.nrows(), A.ncols(), A.get(), A.lda(), contiguous_cache_blocks);
	  seq_outputs[P_first] = results.first;
	  A_top = A;
	}
      else
	{
	  // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
	  const size_t P_mid = (P_first + P_last) / 2;
	  split_t A_split = 
	    partitioner_.split (A, P_first, P_mid, P_last,
				contiguous_cache_blocks);
	  A_top = factor_helper (P_first, P_mid, depth+1, A_split.first, 
				 seq_outputs, par_outputs, R, ldr, 
				 contiguous_cache_blocks);
	  mat_view A_bot = 
	    factor_helper (P_mid+1, P_last, depth+1, A_split.second,
			   seq_outputs, par_outputs, R, ldr, 
			   contiguous_cache_blocks);
	  // Combine the two results
	  factor_pair (P_first, P_mid+1, A_top, A_bot, par_outputs, 
		       contiguous_cache_blocks);
	}

      // If we're completely done, extract the final R factor from
      // the topmost partition.
      if (depth == 0)
	{
#ifdef TBB_DEBUG
	  cerr << "factor_helper: On P_first = " << P_first 
	       << ", extracting R:" << endl
	       << "A_top is " << A_top.nrows() << " x " << A_top.ncols() 
	       << " with leading dimension " << A_top.lda();
#endif // TBB_DEBUG
	  seq_.extract_R (A_top.nrows(), A_top.ncols(), A_top.get(), 
			  A_top.lda(), R, ldr, contiguous_cache_blocks);
	}
      return A_top;
    }


    template< class LocalOrdinal, class Scalar >
    bool
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    apply_helper_empty (const size_t P_first,
			const size_t P_last,
			ConstMatView< LocalOrdinal, Scalar >& Q,
			MatView< LocalOrdinal, Scalar >& C) const
    {
      if (Q.empty())
	{
	  if (! C.empty())
	    throw std::logic_error("Q is empty but C is not!");
	  else
	    return true;
	}
      else if (C.empty())
	{
	  if (! Q.empty())
	    throw std::logic_error("C is empty but Q is not!");
	  else
	    return true;
	}
      else if (P_first > P_last)
	return true;
      else
	return false;
    }


    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    build_partition_array (const size_t P_first,
			   const size_t P_last,
			   typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::array_top_blocks_t& top_blocks,
			   ConstMatView< LocalOrdinal, Scalar >& Q,
			   MatView< LocalOrdinal, Scalar >& C,
			   const bool contiguous_cache_blocks) const
    {
#ifdef TBB_DEBUG
      cerr << "build_partition_array: [" << P_first << ", " << P_last << "]:" << endl
	   << "Q is " << Q.nrows() << " x " << Q.ncols() << " w/ LDA = " 
	   << Q.lda() << endl << "C is " << C.nrows() << " x " << C.ncols() 
	   << " w/ LDA = " << C.lda() << endl;
#endif // TBB_DEBUG

      if (P_first > P_last)
	return;
      else if (P_first == P_last)
	{
	  CacheBlocker< LocalOrdinal, Scalar > blocker (Q.nrows(), Q.ncols(), seq_.cache_blocking_strategy());
	  const_mat_view Q_top = blocker.top_block (Q, contiguous_cache_blocks);
	  mat_view C_top = blocker.top_block (C, contiguous_cache_blocks);
	  top_blocks[P_first] = 
	    std::make_pair (const_mat_view (Q_top.ncols(), Q_top.ncols(), Q_top.get(), Q_top.lda()), 
			    mat_view (C_top.ncols(), C_top.ncols(), C_top.get(), C_top.lda()));
	}
      else
	{
	  // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
	  const size_t P_mid = (P_first + P_last) / 2;
	  const_split_t Q_split = 
	    partitioner_.split (Q, P_first, P_mid, P_last,
				contiguous_cache_blocks);
	  split_t C_split = 
	    partitioner_.split (C, P_first, P_mid, P_last,
				contiguous_cache_blocks);
	  build_partition_array (P_first, P_mid, top_blocks, Q_split.first, 
				 C_split.first, contiguous_cache_blocks);
	  build_partition_array (P_mid+1, P_last, top_blocks, Q_split.second, 
				 C_split.second, contiguous_cache_blocks);
	}
    }
      

    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    apply_helper (const size_t P_first, 
		  const size_t P_last,
		  ConstMatView< LocalOrdinal, Scalar > Q,
		  MatView< LocalOrdinal, Scalar > C,
		  typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::array_top_blocks_t& top_blocks, 
		  const FactorOutput& factor_output,
		  const bool contiguous_cache_blocks) const
    {
      typedef std::pair< const_mat_view, mat_view > apply_t;
#ifdef TBB_DEBUG
      cerr << "apply_helper: [" << P_first << ", " << P_last << "]:" << endl
	   << "Q is " << Q.nrows() << " x " << Q.ncols() << " w/ LDA = " 
	   << Q.lda() << endl << "C is " << C.nrows() << " x " << C.ncols() 
	   << " w/ LDA = " << C.lda() << endl;
#endif // TBB_DEBUG

      if (apply_helper_empty (P_first, P_last, Q, C))
	return;
      else if (P_first == P_last)
	{
	  const std::vector< SeqOutput >& seq_outputs = factor_output.first;
	  seq_.apply ("N", Q.nrows(), Q.ncols(), Q.get(), Q.lda(), 
		      seq_outputs[P_first], C.ncols(), C.get(), 
		      C.lda(), contiguous_cache_blocks);
#ifdef TBB_DEBUG
	  cerr << "BOO!!!" << endl;
#endif // TBB_DEBUG
	}
      else
	{
	  // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
	  const size_t P_mid = (P_first + P_last) / 2;
	  const_split_t Q_split = 
	    partitioner_.split (Q, P_first, P_mid, P_last,
				contiguous_cache_blocks);
	  split_t C_split = 
	    partitioner_.split (C, P_first, P_mid, P_last,
				contiguous_cache_blocks);
	  const ParOutput& par_output = factor_output.second;

	  apply_pair ("N", P_first, P_mid+1, top_blocks[P_mid+1].first,
		      par_output, top_blocks[P_first].second, 
		      top_blocks[P_mid+1].second, contiguous_cache_blocks);
	  apply_helper (P_first, P_mid, Q_split.first, C_split.first,
			top_blocks, factor_output, contiguous_cache_blocks);
	  apply_helper (P_mid+1, P_last, Q_split.second, C_split.second,
			top_blocks, factor_output, contiguous_cache_blocks);
	}
    }


    template< class LocalOrdinal, class Scalar >
    typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::top_blocks_t
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    apply_transpose_helper (const std::string& op,
			    const size_t P_first, 
			    const size_t P_last,
			    ConstMatView< LocalOrdinal, Scalar > Q,
			    MatView< LocalOrdinal, Scalar > C,
			    const typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::FactorOutput& factor_output,
			    const bool contiguous_cache_blocks) const
    {
      if (apply_helper_empty (P_first, P_last, Q, C))
	return std::make_pair (Q, C);
      else if (P_first == P_last)
	{
	  const std::vector< SeqOutput >& seq_outputs = factor_output.first;
	  seq_.apply (op, Q.nrows(), Q.ncols(), Q.get(), Q.lda(), 
		      seq_outputs[P_first], C.ncols(), C.get(), 
		      C.lda(), contiguous_cache_blocks);
	  return std::make_pair (Q, C);
	}
      else
	{
	  // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
	  const size_t P_mid = (P_first + P_last) / 2;

	  const_split_t Q_split = 
	    partitioner_.split (Q, P_first, P_mid, P_last,
				contiguous_cache_blocks);
	  split_t C_split = 
	    partitioner_.split (C, P_first, P_mid, P_last,
				contiguous_cache_blocks);
	  const ParOutput& par_output = factor_output.second;
	  top_blocks_t Top = 
	    apply_transpose_helper (op, P_first, P_mid, Q_split.first, 
				    C_split.first, factor_output, 
				    contiguous_cache_blocks);
	  top_blocks_t Bottom = 
	    apply_transpose_helper (op, P_mid+1, P_last, Q_split.second, 
				    C_split.second, factor_output, 
				    contiguous_cache_blocks);
	  apply_pair (op, P_first, P_mid+1, Bottom.first,
		      par_output, Top.second, Bottom.second,
		      contiguous_cache_blocks);
	  return Top;
	}
    }


    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    factor_pair (const size_t P_top,
		 const size_t P_bot,
		 mat_view& A_top,
		 mat_view& A_bot,
		 std::vector< std::vector< Scalar > >& par_outputs,
		 const bool contiguous_cache_blocks) const
    {
      if (P_top == P_bot) 
	{
	  throw std::logic_error("factor_pair: should never get here!");
	  return; // to pacify the compiler
	}
      // We only read and write the upper ncols x ncols triangle of
      // each block.
      const LocalOrdinal ncols = A_top.ncols();
      if (A_bot.ncols() != ncols)
	throw std::logic_error("A_bot.ncols() != A_top.ncols()");

      std::vector< Scalar >& tau = par_outputs[P_bot];
      std::vector< Scalar > work (ncols);

      TSQR::Combine< LocalOrdinal, Scalar > combine_;
      combine_.factor_pair (ncols, A_top.get(), A_top.lda(),
			    A_bot.get(), A_bot.lda(), &tau[0], &work[0]);
    }

    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    apply_pair (const std::string& trans,
		const size_t P_top,
		const size_t P_bot,
		ConstMatView< LocalOrdinal, Scalar >& Q_bot,
		const std::vector< std::vector< Scalar > >& tau_arrays,
		MatView< LocalOrdinal, Scalar >& C_top,
		MatView< LocalOrdinal, Scalar >& C_bot,
		const bool contiguous_cache_blocks) const
    {
      if (P_top == P_bot) 
	{
	  throw std::logic_error("apply_pair: should never get here!");
	  return; // to pacify the compiler
	}
      const std::vector< Scalar >& tau = tau_arrays[P_bot];
      std::vector< Scalar > work (C_top.ncols());

      TSQR::Combine< LocalOrdinal, Scalar > combine_;
      combine_.apply_pair (trans.c_str(), C_top.ncols(), Q_bot.ncols(), 
			   Q_bot.get(), Q_bot.lda(), &tau[0],
			   C_top.get(), C_top.lda(), 
			   C_bot.get(), C_bot.lda(), &work[0]);
    }

    template< class LocalOrdinal, class Scalar >
    void 
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    cache_block_helper (MatView< LocalOrdinal, Scalar >& A_out,
			ConstMatView< LocalOrdinal, Scalar >& A_in,
			const size_t P_first,
			const size_t P_last) const
    {
      if (P_first > P_last) 
	return;
      else if (P_first == P_last)
	seq_.cache_block (A_out.nrows(), A_out.ncols(), A_out.get(), 
			  A_in.get(), A_in.lda());
      else
	{
	  const size_t P_mid = (P_first + P_last) / 2;
	  const_split_t A_in_split = 
	    partitioner_.split (A_in, P_first, P_mid, P_last, false);
	  split_t A_out_split = 
	    partitioner_.split (A_out, P_first, P_mid, P_last, true);
	  cache_block_helper (A_out_split.first, A_in_split.first, 
			      P_first, P_mid);
	  cache_block_helper (A_out_split.second, A_in_split.second, 
			      P_mid+1, P_last);
	}
    }

    template< class LocalOrdinal, class Scalar >
    void 
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    un_cache_block_helper (MatView< LocalOrdinal, Scalar >& A_out,
			   const ConstMatView< LocalOrdinal, Scalar >& A_in,
			   const size_t P_first,
			   const size_t P_last) const
    {
      if (P_first > P_last) 
	return;
      else if (P_first == P_last)
	seq_.un_cache_block (A_out.nrows(), A_out.ncols(), A_out.get(), 
			     A_out.lda(), A_in.get());
      else
	{
	  const size_t P_mid = (P_first + P_last) / 2;
	  const const_split_t A_in_split = 
	    partitioner_.split (A_in, P_first, P_mid, P_last, true);
	  split_t A_out_split = 
	    partitioner_.split (A_out, P_first, P_mid, P_last, false);
	  
	  un_cache_block_helper (A_out_split.first, A_in_split.first, 
				 P_first, P_mid);
	  un_cache_block_helper (A_out_split.second, A_in_split.second, 
				 P_mid+1, P_last);
	}
    }

    template< class LocalOrdinal, class Scalar >
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    TbbRecursiveTsqr (const size_t num_cores,
		      const size_t cache_size_hint)
      : seq_ (cache_size_hint), ncores_ (1)
    {
      if (num_cores < 1)
	ncores_ = 1; // default is no parallelism
      else
	ncores_ = num_cores;
    }

    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    cache_block (const LocalOrdinal nrows,
		 const LocalOrdinal ncols, 
		 Scalar A_out[],
		 const Scalar A_in[],
		 const LocalOrdinal lda_in) const
    {
      const_mat_view A_in_view (nrows, ncols, A_in, lda_in);
      // Leading dimension doesn't matter, since we're going to cache block it.
      mat_view A_out_view (nrows, ncols, A_out, lda_in);
      cache_block_helper (A_out_view, A_in_view, 0, ncores()-1);
    }

    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    un_cache_block (const LocalOrdinal nrows,
		    const LocalOrdinal ncols,
		    Scalar A_out[],
		    const LocalOrdinal lda_out,		    
		    const Scalar A_in[]) const
    {
      // Leading dimension doesn't matter, since it's cache-blocked.
      const_mat_view A_in_view (nrows, ncols, A_in, lda_out);
      mat_view A_out_view (nrows, ncols, A_out, lda_out);
      un_cache_block_helper (A_out_view, A_in_view, 0, ncores()-1);
    }

    template< class LocalOrdinal, class Scalar >
    typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::FactorOutput
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols, 
	    Scalar A[],
	    const LocalOrdinal lda,
	    Scalar R[],
	    const LocalOrdinal ldr,
	    const bool contiguous_cache_blocks) const
    {
      mat_view A_view (nrows, ncols, A, lda);
      std::vector< SeqOutput > seq_outputs (ncores());
      ParOutput par_outputs (ncores(), std::vector< Scalar >(ncols));
      (void) factor_helper (0, ncores()-1, 0, A_view, seq_outputs, 
			    par_outputs, R, ldr, contiguous_cache_blocks);
      return std::make_pair (seq_outputs, par_outputs);
    }

    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    apply (const std::string& op,
	   const LocalOrdinal nrows,
	   const LocalOrdinal ncols_C,
	   Scalar C[],
	   const LocalOrdinal ldc,
	   const LocalOrdinal ncols_Q,
	   const Scalar Q[],
	   const LocalOrdinal ldq,
	   const typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::FactorOutput& factor_output,
	   const bool contiguous_cache_blocks) const
    {
      const ApplyType apply_type (op);
      if (apply_type == ApplyType::ConjugateTranspose && 
	  ScalarTraits< Scalar >::is_complex)
	throw std::logic_error("Applying Q^H for complex scalar types "
			       "not yet implemented");

      const_mat_view Q_view (nrows, ncols_Q, Q, ldq);
      mat_view C_view (nrows, ncols_C, C, ldc);
      if (! apply_type.transposed())
	{
	  array_top_blocks_t top_blocks (ncores());
	  build_partition_array (0, ncores()-1, top_blocks, Q_view, 
				 C_view, contiguous_cache_blocks);
	  apply_helper (0, ncores()-1, Q_view, C_view, top_blocks, 
			factor_output, contiguous_cache_blocks);
	}
      else
	apply_transpose_helper (op, 0, ncores()-1, Q_view, C_view, 
				factor_output, contiguous_cache_blocks);
    }


    template< class LocalOrdinal, class Scalar >
    void
    TbbRecursiveTsqr< LocalOrdinal, Scalar >::
    explicit_Q (const LocalOrdinal nrows,
		const LocalOrdinal ncols_Q_in,
		const Scalar Q_in[],
		const LocalOrdinal ldq_in,
		const LocalOrdinal ncols_Q_out,
		Scalar Q_out[],
		const LocalOrdinal ldq_out,
		const typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::FactorOutput& factor_output,
		const bool contiguous_cache_blocks) const
    {
      if (ncols_Q_out != ncols_Q_in)
	throw std::logic_error("FIXME Currently, explicit_Q() only works for ncols_Q_out == ncols_Q_in");

      const_mat_view Q_in_view (nrows, ncols_Q_in, Q_in, ldq_in);
      mat_view Q_out_view (nrows, ncols_Q_out, Q_out, ldq_out);

      explicit_Q_helper (0, ncores()-1, Q_out_view, contiguous_cache_blocks);
      apply ("N", nrows, ncols_Q_out, Q_out, ldq_out, ncols_Q_in, 
	     Q_in, ldq_in, factor_output, contiguous_cache_blocks);
    }

  } // namespace TBB
} // namespace TSQR


#endif // __TSQR_TBB_TbbRecursiveTsqr_Def_hpp
