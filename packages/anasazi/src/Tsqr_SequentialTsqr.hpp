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

#ifndef __TSQR_Tsqr_SequentialTsqr_hpp
#define __TSQR_Tsqr_SequentialTsqr_hpp

#include <Tsqr_MatView.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_CacheBlockingStrategy.hpp>
#include <Tsqr_CacheBlocker.hpp>
#include <Tsqr_LocalVerify.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_Combine.hpp>
#include <Tsqr_Util.hpp>

#include <string>
#include <utility> // std::pair
#include <vector>


// #define TSQR_SEQ_TSQR_EXTRA_DEBUG 1

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
#  include <iostream>

template< class MatrixView >
void view_print (const MatrixView& view) {
  using std::cerr;
  using std::endl;
  cerr << view.nrows() << ", " << view.ncols() << ", " << view.lda();
}
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template< class LocalOrdinal, class Scalar >
  class SequentialTsqr {
  private:
    typedef typename std::vector< std::vector< Scalar > >::const_iterator FactorOutputIter;
    typedef typename std::vector< std::vector< Scalar > >::const_reverse_iterator FactorOutputReverseIter;
    typedef MatView< LocalOrdinal, Scalar > mat_view;
    typedef ConstMatView< LocalOrdinal, Scalar > const_mat_view;
    typedef std::pair< mat_view, mat_view > block_pair_type;
    typedef std::pair< const_mat_view, const_mat_view > const_block_pair_type;

    /// Compute the QR factorization in place of the first cache block
    /// A_top in place.  Overwrite the upper triangle of A_top with the
    /// R factor, and return a view of the R factor (stored in place in
    /// A_top).  Overwrite the (strict) lower triangle of A_top, and
    /// the A_top.ncols() entries of tau, with an implicit
    /// representation of the Q factor.  
    ///
    /// \param work [out] Workspace array of length >= A_top.ncols()
    mat_view
    factor_first_block (LAPACK< LocalOrdinal, Scalar >& lapack,
			mat_view& A_top,
			Scalar tau[],
			Scalar work[])
    {
      const LocalOrdinal ncols = A_top.ncols();

      // info must be an int, not a LocalOrdinal, since LAPACK
      // routines always (???) use int for the INFO output argument.
      int info = 0;
      lapack.GEQR2 (A_top.nrows(), A_top.ncols(), A_top.get(), A_top.lda(), tau, work, &info);
      if (info != 0)
	{
	  std::ostringstream os;
	  os << "GEQR2 failed with INFO == " << info;
	  throw std::logic_error (os.str());
	}

      return mat_view(ncols, ncols, A_top.get(), A_top.lda());
    }

    /// Apply the Q factor of the first (topmost) cache blocks, as
    /// computed by factor_first_block() and stored implicitly in
    /// Q_first and tau, to the first (topmost) block C_first of the
    /// matrix C.
    void 
    apply_first_block (const std::string& trans,
		       LAPACK< LocalOrdinal, Scalar >& lapack,
		       const const_mat_view& Q_first,
		       const std::vector< Scalar >& tau,
		       mat_view& C_first,
		       std::vector< Scalar >& work)
    {
      const LocalOrdinal nrows_local = Q_first.nrows();
      int info = 0;
      lapack.ORM2R ("L", trans.c_str(), nrows_local, 
		    C_first.ncols(), Q_first.ncols(),
		    Q_first.get(), Q_first.lda(), &tau[0], 
		    C_first.get(), C_first.lda(), &work[0], &info);
      if (info != 0)
	{
	  std::ostringstream os;
	  os << "ORM2R failed with INFO == " << info;
	  throw std::logic_error (os.str());
	}
    }

    void 
    combine_apply_transpose (Combine< LocalOrdinal, Scalar >& combine,
			     const const_mat_view& Q_cur,
			     const std::vector< Scalar >& tau,
			     mat_view& C_top,
			     mat_view& C_cur,
			     std::vector< Scalar >& work)
    {
      const LocalOrdinal nrows_local = C_cur.nrows();
      combine.apply_inner ("T", nrows_local, 
			   C_cur.ncols(), Q_cur.ncols(),
			   Q_cur.get(), Q_cur.lda(), &tau[0], 
			   C_top.get(), C_top.lda(),
			   C_cur.get(), C_cur.lda(),
			   &work[0]);
    }

    void
    combine_apply (Combine< LocalOrdinal, Scalar >& combine,
		   const const_mat_view& Q_cur,
		   const std::vector< Scalar >& tau,
		   mat_view& C_top,
		   mat_view& C_cur,
		   std::vector< Scalar >& work)
    {
      const LocalOrdinal nrows_local = Q_cur.nrows();
      const LocalOrdinal ncols_Q = Q_cur.ncols();
      const LocalOrdinal ncols_C = C_cur.ncols();

      combine.apply_inner ("N", nrows_local, ncols_C, ncols_Q, 
			   Q_cur.get(), C_cur.lda(), &tau[0],
			   C_top.get(), C_top.lda(),
			   C_cur.get(), C_cur.lda(), &work[0]);
    }

    void
    combine_factor (Combine< LocalOrdinal, Scalar >& combine,
		    mat_view& R,
		    mat_view& A_cur,
		    Scalar tau[],
		    Scalar work[])
    {
      const LocalOrdinal nrows_local = A_cur.nrows();
      const LocalOrdinal ncols = A_cur.ncols();

      combine.factor_inner (nrows_local, ncols, R.get(), R.lda(), 
			    A_cur.get(), A_cur.lda(), tau, work);
    }

  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;
    typedef std::vector< std::vector< Scalar > > FactorOutput;

    /// Constructor
    ///
    /// \param cacheBlockSize [in] Size in bytes of the cache block
    ///   to use in the sequential TSQR factorization.  If 0, the
    ///   implementation will pick a reasonable size.
    SequentialTsqr (const size_t cacheBlockSize = 0) :
      strategy_ (cacheBlockSize) 
    {}

    /// Whether or not the R factor from the QR factorization has a
    /// nonnegative diagonal.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return LAPACK< LocalOrdinal, Scalar >::QR_produces_R_factor_with_nonnegative_diagonal;
    }

    /// \return Cache block size in bytes
    size_t
    cache_block_size () const { return strategy_.cache_block_size(); }

    /// Compute the QR factorization in place of the nrows by ncols
    /// matrix A, with nrows >= ncols, stored either in column-major
    /// order (the default) or as contiguous column-major cache
    /// blocks, with leading dimension lda >= nrows.  Write the
    /// resulting R factor to the top block of A (in place).  (You can
    /// get a view of this via the top_block() method.)  Everything
    /// below the upper triangle of A is overwritten with part of the
    /// implicit representation of the Q factor.  The other part of
    /// that representation is returned.
    FactorOutput
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols,
	    Scalar A[],
	    const LocalOrdinal lda, 
	    const bool contiguous_cache_blocks = false)
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      LAPACK< LocalOrdinal, Scalar > lapack;
      Combine< LocalOrdinal, Scalar > combine;
      std::vector< Scalar > work (ncols);
      FactorOutput tau_arrays;

      // We say "A_rest" because it points to the remaining part of
      // the matrix left to factor; at the beginning, the "remaining"
      // part is the whole matrix, but that will change as the
      // algorithm progresses.
      //
      // Note: if the cache blocks are stored contiguously, lda won't
      // be the correct leading dimension of A, but it won't matter:
      // we only ever operate on A_cur here, and A_cur's leading
      // dimension is set correctly by A_rest.split_top().
      mat_view A_rest (nrows, ncols, A, lda);
      // This call modifies A_rest.
      mat_view A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);

      // Factor the topmost block of A.
      std::vector< Scalar > tau_first (ncols);
      mat_view R_view = factor_first_block (lapack, A_cur, &tau_first[0], &work[0]);
      tau_arrays.push_back (tau_first);

      while (! A_rest.empty())
	{
	  A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	  std::vector< Scalar > tau (ncols);
	  combine_factor (combine, R_view, A_cur, &tau[0], &work[0]);
	  tau_arrays.push_back (tau);
	}
      return tau_arrays;
    }

    /// Extract R factor from factor() results.
    void
    extract_R (const LocalOrdinal nrows,
	       const LocalOrdinal ncols,
	       const Scalar A[],
	       const LocalOrdinal lda,
	       Scalar R[],
	       const LocalOrdinal ldr,
	       const bool contiguous_cache_blocks = false)
    {
      const_mat_view A_view (nrows, ncols, A, lda);

      // Identify top cache block of A
      const_mat_view A_top = top_block (A_view, contiguous_cache_blocks);

      // Fill R (including lower triangle) with zeros.
      fill_matrix (ncols, ncols, R, ldr, Scalar(0));

      // Copy out the upper triangle of the R factor from A into R.
      copy_upper_triangle (ncols, ncols, R, ldr, A_top.get(), A_top.lda());
    }

    /// Compute the QR factorization in place of the nrows by ncols
    /// matrix A, with nrows >= ncols, stored either in column-major
    /// order (the default) or as contiguous column-major cache
    /// blocks, with leading dimension lda >= nrows.  Return an
    /// implicit representation of the Q factor.  Copy the resulting R
    /// factor to the R output argument.
    ///
    /// \note An alternate syntax for factor().  Useful when using
    /// SequentialTsqr by itself, rather than (say) in TbbTsqr.
    FactorOutput
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols,
	    Scalar A[],
	    const LocalOrdinal lda, 
	    Scalar R[],
	    const LocalOrdinal ldr,
	    const bool contiguous_cache_blocks = false)
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      LAPACK< LocalOrdinal, Scalar > lapack;
      Combine< LocalOrdinal, Scalar > combine;
      std::vector< Scalar > work (ncols);
      FactorOutput tau_arrays;

      // We say "A_rest" because it points to the remaining part of
      // the matrix left to factor; at the beginning, the "remaining"
      // part is the whole matrix, but that will change as the
      // algorithm progresses.
      //
      // Note: if the cache blocks are stored contiguously, lda won't
      // be the correct leading dimension of A, but it won't matter:
      // we only ever operate on A_cur here, and A_cur's leading
      // dimension is set correctly by A_rest.split_top().
      mat_view A_rest (nrows, ncols, A, lda);
      // This call modifies A_rest.
      mat_view A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);

      // Factor the topmost block of A.
      std::vector< Scalar > tau_first (ncols);
      mat_view R_view = factor_first_block (lapack, A_cur, &tau_first[0], &work[0]);
      tau_arrays.push_back (tau_first);

      while (! A_rest.empty())
	{
	  A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	  std::vector< Scalar > tau (ncols);
	  combine_factor (combine, R_view, A_cur, &tau[0], &work[0]);
	  tau_arrays.push_back (tau);
	}
      
      // Copy the R factor resulting from the factorization out of
      // R_view (a view of the topmost cache block of A) into the R
      // output argument.
      fill_matrix (ncols, ncols, R, ldr, Scalar(0));
      copy_upper_triangle (ncols, ncols, R, ldr, R_view.get(), R_view.lda());
      return tau_arrays;
    }

    void
    apply (const std::string& op,
	   const LocalOrdinal nrows,
	   const LocalOrdinal ncols_Q,
	   const Scalar* const Q,
	   const LocalOrdinal ldq,
	   const FactorOutput& factor_output,
	   const LocalOrdinal ncols_C,
	   Scalar* const C,
	   const LocalOrdinal ldc,
	   const bool contiguous_cache_blocks = false)
    {
      // Quick exit and error tests
      if (ncols_Q == 0 || ncols_C == 0 || nrows == 0)
	return;
      else if (ldc < nrows)
	{
	  std::ostringstream os;
	  os << "SequentialTsqr::apply: ldc (= " << ldc << ") < nrows (= " << nrows << ")";
	  throw std::invalid_argument (os.str());
	}
      else if (ldq < nrows)
	{
	  std::ostringstream os;
	  os << "SequentialTsqr::apply: ldq (= " << ldq << ") < nrows (= " << nrows << ")";
	  throw std::invalid_argument (os.str());
	}

      // If contiguous cache blocks are used, then we have to use the
      // same convention as we did for factor().  Otherwise, we are
      // free to choose the cache block dimensions as we wish in
      // apply(), independently of what we did in factor().
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols_Q, strategy_);
      LAPACK< LocalOrdinal, Scalar > lapack;
      Combine< LocalOrdinal, Scalar > combine;

      bool transposed;
      if (op[0] == 'H' || op[0] == 'h')
	throw std::logic_error("TSQR::apply: applying Q^H not yet implemented");
      else if (op[0] == 'N' || op[0] == 'n')
	transposed = false;
      else if (op[0] == 't' || op[0] == 't')
	transposed = true;
      else
	throw std::invalid_argument ("SequentialTsqr::apply: Invalid op argument \"" + op + "\"");

      const FactorOutput& tau_arrays = factor_output; // rename for encapsulation
      std::vector< Scalar > work (ncols_C);
      
      // We say "*_rest" because it points to the remaining part of
      // the matrix left to factor; at the beginning, the "remaining"
      // part is the whole matrix, but that will change as the
      // algorithm progresses.
      //
      // Note: if the cache blocks are stored contiguously, ldq
      // resp. ldc won't be the correct leading dimension, but it
      // won't matter, since we only read the leading dimension of
      // return values of split_top_block() / split_bottom_block(),
      // which are set correctly (based e.g., on whether or not we are
      // using contiguous cache blocks).
      const_mat_view Q_rest (nrows, ncols_Q, Q, ldq);
      mat_view C_rest (nrows, ncols_C, C, ldc);

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
      if (Q_rest.empty())
	throw std::logic_error ("Q_rest initially empty");
      else if (C_rest.empty())
	throw std::logic_error ("C_rest initially empty");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG

      // Identify the top ncols_C by ncols_C block of C.  C_rest is
      // not modified.
      mat_view C_top = blocker.top_block (C_rest, contiguous_cache_blocks);

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
      if (C_top.empty())
	throw std::logic_error ("C_top initially empty");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG

      if (transposed)
	{
	  const_mat_view Q_cur = blocker.split_top_block (Q_rest, contiguous_cache_blocks);
	  mat_view C_cur = blocker.split_top_block (C_rest, contiguous_cache_blocks);

	  // Apply the topmost block of Q.
	  FactorOutputIter tau_iter = tau_arrays.begin();
	  const std::vector< Scalar >& tau = *tau_iter++;
	  apply_first_block ("T", lapack, Q_cur, tau, C_cur, work);

	  while (! Q_rest.empty())
	    {
	      Q_cur = blocker.split_top_block (Q_rest, contiguous_cache_blocks);
	      C_cur = blocker.split_top_block (C_rest, contiguous_cache_blocks);
#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
	      if (tau_iter == tau_arrays.end())
		throw std::logic_error("Not enough tau arrays!");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG
	      combine_apply_transpose (combine, Q_cur, *tau_iter++, C_top, C_cur, work);
	    }
#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
	  if (tau_iter != tau_arrays.end())
	    throw std::logic_error ("Too many tau arrays!");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG
	}
      else
	{
	  // Start with the last local Q factor and work backwards up the matrix.
	  FactorOutputReverseIter tau_iter = tau_arrays.rbegin();

	  const_mat_view Q_cur = blocker.split_bottom_block (Q_rest, contiguous_cache_blocks);
	  mat_view C_cur = blocker.split_bottom_block (C_rest, contiguous_cache_blocks);

	  while (! Q_rest.empty())
	    {
#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
	      if (Q_cur.empty())
	      	throw std::logic_error ("Q_cur empty at last stage of applying Q");
	      else if (C_cur.empty())
	      	throw std::logic_error ("C_cur empty at last stage of applying Q");
	      else if (tau_iter == tau_arrays.rend())
	      	throw std::logic_error ("Not enough tau arrays!");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG
	      combine_apply (combine, Q_cur, *tau_iter++, C_top, C_cur, work);
	      Q_cur = blocker.split_bottom_block (Q_rest, contiguous_cache_blocks);
	      C_cur = blocker.split_bottom_block (C_rest, contiguous_cache_blocks);
	    }
	  //
	  // Apply to last (topmost) cache block.
	  //
#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
	  if (Q_cur.empty())
	    throw std::logic_error ("Q_cur empty at last stage of applying Q");
	  else if (C_cur.empty())
	    throw std::logic_error ("C_cur empty at last stage of applying Q");
	  else if (tau_iter == tau_arrays.rend())
	    throw std::logic_error ("Not enough tau arrays!");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG

	  apply_first_block ("N", lapack, Q_cur, *tau_iter++, C_cur, work);

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
	  if (tau_iter != tau_arrays.rend())
	    throw std::logic_error ("Too many tau arrays!");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG
	}
    }

    void
    explicit_Q (const LocalOrdinal nrows,
		const LocalOrdinal ncols_Q,
		const Scalar Q[],
		const LocalOrdinal ldq,
		const FactorOutput& factor_output,
		const LocalOrdinal ncols_C,
		Scalar C[],
		const LocalOrdinal ldc,
		const bool contiguous_cache_blocks = false)
    {
      // Identify top ncols_C by ncols_C block of C.  C_view is not
      // modified.  top_block() will set C_top to have the correct
      // leading dimension, whether or not cache blocks are stored
      // contiguously.
      mat_view C_view (nrows, ncols_C, C, ldc);
      mat_view C_top = top_block (C_view, contiguous_cache_blocks);

      // Fill C with zeros, and then fill the topmost block of C with
      // the first ncols_C columns of the identity matrix, so that C
      // itself contains the first ncols_C columns of the identity
      // matrix.
      fill_with_zeros (nrows, ncols_C, C, ldc, contiguous_cache_blocks);
      for (LocalOrdinal j = 0; j < ncols_C; j++)
      	C_top(j, j) = Scalar(1);

      // Apply the Q factor to C, to extract the first ncols_C columns
      // of Q in explicit form.
      apply ("N", nrows, ncols_Q, Q, ldq, factor_output, 
	     ncols_C, C, ldc, contiguous_cache_blocks);
    }


    /// Cache-block the given A_in matrix, writing the results to A_out.
    void
    cache_block (const LocalOrdinal nrows,
		 const LocalOrdinal ncols,
		 Scalar A_out[],
		 const Scalar A_in[],
		 const LocalOrdinal lda_in) const
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      blocker.cache_block (nrows, ncols, A_out, A_in, lda_in);
    }

    /// "Un"-cache-block the given A_in matrix, writing the results to A_out.
    void
    un_cache_block (const LocalOrdinal nrows,
		    const LocalOrdinal ncols,
		    Scalar A_out[],
		    const LocalOrdinal lda_out,		    
		    const Scalar A_in[]) const
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      blocker.un_cache_block (nrows, ncols, A_out, lda_out, A_in);
    }

    /// Fill the nrows by ncols matrix A with zeros.
    void
    fill_with_zeros (const LocalOrdinal nrows,
		     const LocalOrdinal ncols,
		     Scalar A[],
		     const LocalOrdinal lda, 
		     const bool contiguous_cache_blocks = false)
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      blocker.fill_with_zeros (nrows, ncols, A, lda, contiguous_cache_blocks);
    }

    /// Return a view of the topmost cache block (on this node) of the
    /// given matrix C.  NOTE that this is not necessarily square,
    /// though it must have at least as many rows as columns.  For a
    /// square ncols by ncols block, as needed in TSQR::Tsqr::apply(),
    /// if the output is ret, do MatView< LocalOrdinal, Scalar >
    /// (ncols, ncols, ret.get(), ret.lda()) to get an ncols by ncols
    /// block.
    template< class MatrixViewType >
    MatrixViewType
    top_block (const MatrixViewType& C, 
	       const bool contiguous_cache_blocks = false) const 
    {
      // The CacheBlocker object knows how to construct a view of the
      // top cache block of C.  This is complicated because cache
      // blocks (in C) may or may not be stored contiguously.  If they
      // are stored contiguously, the CacheBlocker knows the right
      // layout, based on the cache blocking strategy.
      CacheBlocker< LocalOrdinal, Scalar > blocker (C.nrows(), C.ncols(), strategy_);

      // C_top_block is a view of the topmost cache block of C.
      // C_top_block should have >= ncols rows, otherwise either cache
      // blocking is broken or the input matrix C itself had fewer
      // rows than columns.
      MatrixViewType C_top_block = blocker.top_block (C, contiguous_cache_blocks);
      if (C_top_block.nrows() < C_top_block.ncols())
	throw std::logic_error ("C\'s topmost cache block has fewer rows than "
				"columns");
      return C_top_block;
    }

  private:
    CacheBlockingStrategy< LocalOrdinal, Scalar > strategy_;
  };
  
} // namespace TSQR

#endif // __TSQR_Tsqr_SequentialTsqr_hpp
