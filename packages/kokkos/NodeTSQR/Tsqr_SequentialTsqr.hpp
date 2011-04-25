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

#ifndef __TSQR_Tsqr_SequentialTsqr_hpp
#define __TSQR_Tsqr_SequentialTsqr_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_CacheBlockingStrategy.hpp>
#include <Tsqr_CacheBlocker.hpp>
#include <Tsqr_LocalVerify.hpp>
#include <Tsqr_Combine.hpp>
#include <Tsqr_Util.hpp>

#include <Teuchos_Describable.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <utility> // std::pair
#include <vector>

// #define TSQR_SEQ_TSQR_EXTRA_DEBUG 1

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
#  include <iostream>
using std::cerr;
using std::endl;

template< class MatrixView >
void view_print (const MatrixView& view) {
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
    factor_first_block (Combine< LocalOrdinal, Scalar >& combine,
			mat_view& A_top,
			std::vector< Scalar >& tau,
			std::vector< Scalar >& work)
    {
      const LocalOrdinal ncols = A_top.ncols();
      combine.factor_first (A_top.nrows(), ncols, A_top.get(), A_top.lda(), 
			    &tau[0], &work[0]);
      return mat_view(ncols, ncols, A_top.get(), A_top.lda());
    }

    /// Apply the Q factor of the first (topmost) cache blocks, as
    /// computed by factor_first_block() and stored implicitly in
    /// Q_first and tau, to the first (topmost) block C_first of the
    /// matrix C.
    void 
    apply_first_block (Combine< LocalOrdinal, Scalar >& combine,
		       const ApplyType& applyType,
		       const const_mat_view& Q_first,
		       const std::vector< Scalar >& tau,
		       mat_view& C_first,
		       std::vector< Scalar >& work)
    {
      const LocalOrdinal nrowsLocal = Q_first.nrows();
      combine.apply_first (applyType, nrowsLocal, C_first.ncols(), 
			   Q_first.ncols(), Q_first.get(), Q_first.lda(),
			   &tau[0], C_first.get(), C_first.lda(), &work[0]);
    }

    void
    combine_apply (Combine< LocalOrdinal, Scalar >& combine,
		   const ApplyType& apply_type,
		   const const_mat_view& Q_cur,
		   const std::vector< Scalar >& tau,
		   mat_view& C_top,
		   mat_view& C_cur,
		   std::vector< Scalar >& work)
    {
      const LocalOrdinal nrows_local = Q_cur.nrows();
      const LocalOrdinal ncols_Q = Q_cur.ncols();
      const LocalOrdinal ncols_C = C_cur.ncols();

      combine.apply_inner (apply_type, 
			   nrows_local, ncols_C, ncols_Q, 
			   Q_cur.get(), C_cur.lda(), &tau[0],
			   C_top.get(), C_top.lda(),
			   C_cur.get(), C_cur.lda(), &work[0]);
    }

    void
    combine_factor (Combine< LocalOrdinal, Scalar >& combine,
		    mat_view& R,
		    mat_view& A_cur,
		    std::vector< Scalar >& tau,
		    std::vector< Scalar >& work)
    {
      const LocalOrdinal nrows_local = A_cur.nrows();
      const LocalOrdinal ncols = A_cur.ncols();

      combine.factor_inner (nrows_local, ncols, R.get(), R.lda(), 
			    A_cur.get(), A_cur.lda(), &tau[0], 
			    &work[0]);
    }

  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
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

    /// \brief One-line description of this object.
    ///
    /// This implements Teuchos::Describable::description().  For now,
    /// SequentialTsqr uses the default implementation of
    /// Teuchos::Describable::describe().
    std::string description () const {
      std::ostringstream os;
      os << "Intranode Tall Skinny QR (TSQR): sequential cache-blocked "
	"implementation with cache block size " << this->cache_block_size() 
	 << " bytes.";
      return os.str();
    }

    /// Whether or not the R factor from the QR factorization has a
    /// nonnegative diagonal.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      typedef Combine< LocalOrdinal, Scalar > combine_type;
      return combine_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    //! Cache block size in bytes.
    size_t cache_block_size () const { 
      return strategy_.cache_block_size(); 
    }

    /// \brief Compute QR factorization (implicitly stored Q factor) of A.
    ///
    /// Compute the QR factorization in place of the nrows by ncols
    /// matrix A, with nrows >= ncols, stored either in column-major
    /// order (the default) or as contiguous column-major cache
    /// blocks, with leading dimension lda >= nrows.  Write the
    /// resulting R factor to the top block of A (in place).  (You can
    /// get a view of this via the top_block() method.)  Everything
    /// below the upper triangle of A is overwritten with part of the
    /// implicit representation of the Q factor.  The other part of
    /// that representation is returned.
    ///
    /// \param nrows [in] Number of rows in A
    /// \param ncols [in] Number of columns in A 
    /// \param A [in/out] On input: nrows by ncols dense matrix to
    ///   factor.  On output: partial representation of the implicitly
    ///   stored Q factor.
    /// \param lda [in] Leading dimension of A, if A is stored in
    ///   column-major order.  Otherwise not read.
    /// \param contiguous_cache_blocks [in] Whether the matrix A is
    ///   stored in a contiguously cache-blocked format.
    ///
    /// \return Partial representation of the implicitly stored Q
    ///   factor.  The complete representation includes A (on output).
    ///   The FactorOutput and A go together.
    FactorOutput
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols,
	    Scalar A[],
	    const LocalOrdinal lda, 
	    const bool contiguous_cache_blocks = false)
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
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
      mat_view R_view = factor_first_block (combine, A_cur, tau_first, work);
      tau_arrays.push_back (tau_first);

      while (! A_rest.empty())
	{
	  A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	  std::vector< Scalar > tau (ncols);
	  combine_factor (combine, R_view, A_cur, tau, work);
	  tau_arrays.push_back (tau);
	}
      return tau_arrays;
    }

    //! Extract R factor from factor() results.
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
      mat_view R_view = factor_first_block (combine, A_cur, tau_first, work);
      tau_arrays.push_back (tau_first);

      while (! A_rest.empty())
	{
	  A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	  std::vector< Scalar > tau (ncols);
	  combine_factor (combine, R_view, A_cur, tau, work);
	  tau_arrays.push_back (tau);
	}
      
      // Copy the R factor resulting from the factorization out of
      // R_view (a view of the topmost cache block of A) into the R
      // output argument.
      fill_matrix (ncols, ncols, R, ldr, Scalar(0));
      copy_upper_triangle (ncols, ncols, R, ldr, R_view.get(), R_view.lda());
      return tau_arrays;
    }


    //! Number of cache blocks that factor() would use.
    LocalOrdinal
    factor_num_cache_blocks (const LocalOrdinal nrows,
			     const LocalOrdinal ncols,
			     Scalar A[],
			     const LocalOrdinal lda, 
			     const bool contiguous_cache_blocks = false)
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      LocalOrdinal count = 0;

      mat_view A_rest (nrows, ncols, A, lda);
      if (A_rest.empty())
	return count;

      mat_view A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
      count++; // first factor step

      while (! A_rest.empty())
	{
	  A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	  count++; // next factor step
	}
      return count;
    }

    //! Apply implicit Q factor stored in Q and factor_output to C.
    void
    apply (const ApplyType& apply_type,
	   const LocalOrdinal nrows,
	   const LocalOrdinal ncols_Q,
	   const Scalar Q[],
	   const LocalOrdinal ldq,
	   const FactorOutput& factor_output,
	   const LocalOrdinal ncols_C,
	   Scalar C[],
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

      const bool transposed = apply_type.transposed();
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
	  apply_first_block (combine, apply_type, Q_cur, tau, C_cur, work);

	  while (! Q_rest.empty())
	    {
	      Q_cur = blocker.split_top_block (Q_rest, contiguous_cache_blocks);
	      C_cur = blocker.split_top_block (C_rest, contiguous_cache_blocks);
#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
	      if (tau_iter == tau_arrays.end())
		throw std::logic_error("Not enough tau arrays!");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG
	      combine_apply (combine, apply_type, Q_cur, *tau_iter++, C_top, C_cur, work);
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
	      combine_apply (combine, apply_type, Q_cur, *tau_iter++, C_top, C_cur, work);
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

	  apply_first_block (combine, apply_type, Q_cur, *tau_iter++, C_cur, work);

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
	  if (tau_iter != tau_arrays.rend())
	    throw std::logic_error ("Too many tau arrays!");
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG
	}
    }

    //! Compute explicit Q factor from result of factor().
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
      for (LocalOrdinal j = 0; j < ncols_C; ++j)
      	C_top(j, j) = Scalar(1);

      // Apply the Q factor to C, to extract the first ncols_C columns
      // of Q in explicit form.
      apply (ApplyType::NoTranspose, 
	     nrows, ncols_Q, Q, ldq, factor_output, 
	     ncols_C, C, ldc, contiguous_cache_blocks);
    }


    /// \brief Compute Q*B.
    ///
    /// Compute matrix-matrix product Q*B, where Q is nrows by ncols
    /// and B is ncols by ncols.  Respect cache blocks of Q.
    void
    Q_times_B (const LocalOrdinal nrows,
	       const LocalOrdinal ncols,
	       Scalar Q[],
	       const LocalOrdinal ldq,
	       const Scalar B[],
	       const LocalOrdinal ldb,
	       const bool contiguous_cache_blocks = false) const
    {
      // We don't do any other error checking here (e.g., matrix
      // dimensions), though it would be a good idea to do so.

      // Take the easy exit if available.
      if (ncols == 0 || nrows == 0)
	return;

      // Compute Q := Q*B by iterating through cache blocks of Q.
      // This iteration works much like iteration through cache blocks
      // of A in factor() (which see).  Here, though, each cache block
      // computation is completely independent of the others; a slight
      // restructuring of this code would parallelize nicely using
      // OpenMP.
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      BLAS< LocalOrdinal, Scalar > blas;
      mat_view Q_rest (nrows, ncols, Q, ldq);
      Matrix< LocalOrdinal, Scalar > 
	Q_cur_copy (LocalOrdinal(0), LocalOrdinal(0)); // will be resized
      while (! Q_rest.empty())
	{
	  mat_view Q_cur = 
	    blocker.split_top_block (Q_rest, contiguous_cache_blocks);

	  // GEMM doesn't like aliased arguments, so we use a copy.
	  // We only copy the current cache block, rather than all of
	  // Q; this saves memory.
	  Q_cur_copy.reshape (Q_cur.nrows(), ncols);
	  Q_cur_copy.copy (Q_cur);
	  // Q_cur := Q_cur_copy * B.
	  blas.GEMM ("N", "N", Q_cur.nrows(), ncols, ncols, Scalar(1),
		     Q_cur_copy.get(), Q_cur_copy.lda(), B, ldb,
		     Scalar(0), Q_cur.get(), Q_cur.lda());
	}
    }

    /// \brief Reveal rank of the R factor from \c factor().
    ///
    /// Compute SVD \f$R = U \Sigma V^*\f$, not in place.  Use the
    /// resulting singular values to compute the numerical rank of R,
    /// with respect to the relative tolerance tol.  If R is full
    /// rank, return without modifying R.  If R is not full rank,
    /// overwrite R with \f$\Sigma \cdot V^*\f$.
    ///
    /// \return Numerical rank of R: 0 <= rank <= ncols.
    LocalOrdinal
    reveal_R_rank (const LocalOrdinal ncols,
		   Scalar R[],
		   const LocalOrdinal ldr,
		   Scalar U[],
		   const LocalOrdinal ldu,
		   const magnitude_type tol) const 
    {
      if (tol < 0)
	{
	  std::ostringstream os;
	  os << "reveal_R_rank: negative tolerance tol = "
	     << tol << " is not allowed.";
	  throw std::logic_error (os.str());
	}
      // We don't do any other error checking here (e.g., matrix
      // dimensions), though it would be a good idea to do so.

      // Take the easy exit if available.
      if (ncols == 0)
	return 0;

      LAPACK< LocalOrdinal, Scalar > lapack;
      MatView< LocalOrdinal, Scalar > R_view (ncols, ncols, R, ldr);
      Matrix< LocalOrdinal, Scalar > B (R_view); // B := R (deep copy)
      MatView< LocalOrdinal, Scalar > U_view (ncols, ncols, U, ldu);
      Matrix< LocalOrdinal, Scalar > VT (ncols, ncols, Scalar(0));

      std::vector< magnitude_type > svd_rwork (5*ncols);
      std::vector< magnitude_type > singular_values (ncols);
      LocalOrdinal svd_lwork = -1; // -1 for LWORK query; will be changed
      int svd_info = 0;

      // LAPACK LWORK query for singular value decomposition.  WORK
      // array is always of ScalarType, even in the complex case.
      {
	Scalar svd_lwork_scalar = Scalar(0);
	lapack.GESVD ("A", "A", ncols, ncols, 
		       B.get(), B.lda(), &singular_values[0], 
		       U_view.get(), U_view.lda(), VT.get(), VT.lda(),
		       &svd_lwork_scalar, svd_lwork, &svd_rwork[0], &svd_info);
	if (svd_info != 0)
	  {
	    std::ostringstream os;
	    os << "reveal_R_rank: GESVD LWORK query returned nonzero INFO = "
	       << svd_info;
	    throw std::logic_error (os.str());
	  }
	// Scalar could be a complex type, but LAPACK should only ever
	// return a Scalar with zero imaginary part.
	if (Teuchos::ScalarTraits< Scalar >::isComplex && Teuchos::ScalarTraits< Scalar >::imag (svd_lwork_scalar) != Teuchos::ScalarTraits< magnitude_type >::zero())
	  {
	    std::ostringstream os;
	    os << "In SequentialTsqr::reveal_rank: GESVD LWORK query returned "
	      "LWORK (as a Scalar) with a nonzero imaginary part.  How do I "
	      "allocate an imaginary amount of workspace?  LAPACK should never "
	      "do this, so there must be some weird internal corruption "
	      "somewhere.  Returned LWORK value is " << svd_lwork_scalar << ".";
	    throw std::logic_error (os.str());
	  }
	// Scalar has a zero imaginary part, so we can try converting
	// its real part into a LocalOrdinal.
	svd_lwork = static_cast< LocalOrdinal > (Teuchos::ScalarTraits< Scalar >::real (svd_lwork_scalar));
	// Check the LWORK cast.  LAPACK shouldn't ever return LWORK
	// that won't fit in an OrdinalType, but it's not bad to make
	// sure.
	if (static_cast< Scalar > (svd_lwork) != svd_lwork_scalar)
	  {
	    std::ostringstream os;
	    os << "In SequentialTsqr::reveal_rank: GESVD LWORK query "
	      "returned LWORK that doesn\'t fit in LocalOrdinal: returned "
	      "LWORK (as Scalar) is " << svd_lwork_scalar << ", but cast to "
	      "LocalOrdinal it becomes " << svd_lwork << ".";
	    throw std::logic_error (os.str());
	  }
	// Make sure svd_lwork >= 0.
	if (std::numeric_limits< LocalOrdinal >::is_signed && svd_lwork < 0)
	  {
	    std::ostringstream os;
	    os << "In SequentialTsqr::reveal_rank: GESVD LWORK query "
	      "returned negative LWORK = " << svd_lwork;
	    throw std::logic_error (os.str());
	  }
      }
      // Allocate workspace for SVD.
      std::vector< Scalar > svd_work (svd_lwork);

      // Compute SVD $B := U \Sigma V^*$.  B is overwritten, which is
      // why we copied R into B (so that we don't overwrite R if R is
      // full rank).
      lapack.GESVD ("A", "A", ncols, ncols, 
		    B.get(), B.lda(), &singular_values[0], 
		    U_view.get(), U_view.lda(), VT.get(), VT.lda(),
		    &svd_work[0], svd_lwork, &svd_rwork[0], &svd_info);

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
      {
      	cerr << "-- GESVD computed singular values:" << endl;
      	for (int k = 0; k < ncols; ++k)
      	  cerr << singular_values[k] << " ";
      	cerr << endl;
      }
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG

      // GESVD computes singular values in decreasing order and they
      // are all nonnegative.  We know by now that ncols > 0.  "tol"
      // is a relative tolerance: relative to the largest singular
      // value, which is the 2-norm of the matrix.
      const magnitude_type absolute_tolerance = tol * singular_values[0];

      // Determine rank of B, using singular values.  
      LocalOrdinal rank = ncols;
      for (LocalOrdinal k = 1; k < ncols; ++k)
	// "<=" in case singular_values[0] == 0.
	if (singular_values[k] <= absolute_tolerance)
	  {
	    rank = k;
	    break;
	  }

      if (rank == ncols)
	return rank; // Don't modify Q or R, if R is full rank.

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
      {
	cerr << "Rank of B (i.e., R): " << rank << " < ncols=" << ncols << endl;
	cerr << "Original R = " << endl;
	print_local_matrix (cerr, ncols, ncols, R, ldr);
      }
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG

      //
      // R is not full rank.  
      //
      // 1. Compute \f$R := \Sigma V^*\f$.
      // 2. Return rank (0 <= rank < ncols).
      //

      // Compute R := \Sigma VT.  \Sigma is diagonal so we apply it
      // column by column (normally one would think of this as row by
      // row, but this "Hadamard product" formulation iterates more
      // efficiently over VT).  
      //
      // After this computation, R may no longer be upper triangular.
      // R may be zero if all the singular values are zero, but we
      // don't need to check for this case; it's rare in practice, and
      // the computations below will be correct regardless.
      for (LocalOrdinal j = 0; j < ncols; ++j)
	{
	  const Scalar* const VT_j = &VT(0,j);
	  Scalar* const R_j = &R_view(0,j);

	  for (LocalOrdinal i = 0; i < ncols; ++i)
	    R_j[i] = singular_values[i] * VT_j[i];
	}

#ifdef TSQR_SEQ_TSQR_EXTRA_DEBUG
      {
	cerr << "Resulting R = " << endl;
	print_local_matrix (cerr, ncols, ncols, R, ldr);
      }
#endif // TSQR_SEQ_TSQR_EXTRA_DEBUG

      return rank;
    }

    /// \brief Rank-revealing decomposition.
    ///
    /// Using the R factor from factor() and the explicit Q factor
    /// from explicit_Q(), compute the SVD of R (\f$R = U \Sigma
    /// V^*\f$).  R.  If R is full rank (with respect to the given
    /// relative tolerance tol), don't change Q or R.  Otherwise,
    /// compute \f$Q := Q \cdot U\f$ and \f$R := \Sigma V^*\f$ in
    /// place (the latter may be no longer upper triangular).
    ///
    /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
    ///
    LocalOrdinal
    reveal_rank (const LocalOrdinal nrows,
		 const LocalOrdinal ncols,
		 Scalar Q[],
		 const LocalOrdinal ldq,
		 Scalar R[],
		 const LocalOrdinal ldr,
		 const magnitude_type tol,
		 const bool contiguous_cache_blocks = false) const
    {
      // Take the easy exit if available.
      if (ncols == 0)
	return 0;
      Matrix< LocalOrdinal, Scalar > U (ncols, ncols, Scalar(0));
      const LocalOrdinal rank = 
	reveal_R_rank (ncols, R, ldr, U.get(), U.ldu(), tol);
      
      if (rank < ncols)
	{
	  // If R is not full rank: reveal_R_rank() already computed
	  // the SVD \f$R = U \Sigma V^*\f$ of (the input) R, and
	  // overwrote R with \f$\Sigma V^*\f$.  Now, we compute \f$Q
	  // := Q \cdot U\f$, respecting cache blocks of Q.
	  Q_times_B (nrows, ncols, Q, ldq, U.get(), U.lda(), 
		     contiguous_cache_blocks);
	}
      return rank;
    }

    /// \brief Cache block A_in into A_out.
    ///
    /// \param nrows [in] Number of rows in A_in and A_out.
    /// \param ncols [in] Number of columns in A_in and A_out.
    /// \param A_out [out] Result of cache-blocking A_in.
    /// \param A_in [in] Matrix to cache block, stored in column-major
    ///   order with leading dimension lda_in.
    /// \param lda_in [in] Leading dimension of A_in.  (See the LAPACK
    ///   documentation for a definition of "leading dimension.")
    ///   lda_in >= nrows.
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

    /// \brief Un - cache block A_in into A_out.
    ///
    /// A_in is a matrix produced by \c cache_block().  It is
    /// organized as contiguously stored cache blocks.  This method
    /// reorganizes A_in into A_out as an ordinary matrix stored in
    /// column-major order with leading dimension lda_out.
    ///
    /// \param nrows [in] Number of rows in A_in and A_out.
    /// \param ncols [in] Number of columns in A_in and A_out.
    /// \param A_out [out] Result of un-cache-blocking A_in.
    ///   Matrix stored in column-major order with leading
    ///   dimension lda_out.
    /// \param lda_out [in] Leading dimension of A_out.  (See the
    ///   LAPACK documentation for a definition of "leading
    ///   dimension.")  lda_out >= nrows.
    /// \param A_in [in] Matrix to un-cache-block.
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

    /// \brief Fill the nrows by ncols matrix A with zeros.
    /// 
    /// Fill the matrix A with zeros, in a way that respects the cache
    /// blocking scheme.
    ///
    /// \param nrows [in] Number of rows in A
    /// \param ncols [in] Number of columns in A 
    /// \param A [out] nrows by ncols column-major-order dense matrix 
    ///   with leading dimension lda
    /// \param lda [in] Leading dimension of A: lda >= nrows
    /// \param contiguous_cache_blocks [in] Whether the cache blocks
    ///   in A are stored contiguously.
    void
    fill_with_zeros (const LocalOrdinal nrows,
		     const LocalOrdinal ncols,
		     Scalar A[],
		     const LocalOrdinal lda, 
		     const bool contiguous_cache_blocks = false) const
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      blocker.fill_with_zeros (nrows, ncols, A, lda, contiguous_cache_blocks);
    }

    /// \brief Return topmost cache block of C
    ///
    /// \param C [in] Matrix (view), supporting the usual nrows(),
    ///   ncols(), get(), lda() interface.
    /// \param contiguous_cache_blocks [in] Whether the cache blocks
    ///   in C are stored contiguously.
    ///
    /// Return a view of the topmost cache block (on this node) of the
    /// given matrix C.  This is not necessarily square, though it
    /// must have at least as many rows as columns.  For a square
    /// ncols by ncols block, as needed by Tsqr::apply(), do as 
    /// follows:
    /// \code 
    /// MatrixViewType top = this->top_block (C, contig);
    /// MatView<LocalOrdinal, Scalar> square (ncols, ncols, top.get(), top.lda());
    /// \endcode
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
      CacheBlocker<LocalOrdinal, Scalar> blocker (C.nrows(), C.ncols(), strategy_);

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
    //! Strategy object that helps us cache block matrices.
    CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
  };
  
} // namespace TSQR

#endif // __TSQR_Tsqr_SequentialTsqr_hpp
