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
#include <Tsqr_Combine.hpp>
#include <Tsqr_LocalVerify.hpp>
#include <Tsqr_NodeTsqr.hpp>
#include <Tsqr_Util.hpp>

#include <Teuchos_Describable.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <utility> // std::pair
#include <vector>

namespace TSQR {

  /// \class SequentialTsqr
  /// \brief Sequential TSQR factorization.
  /// \author Mark Hoemmen
  ///
  /// TSQR (Tall Skinny QR) is a collection of different algorithms
  /// for computing the QR factorization of a "tall and skinny" matrix
  /// (with many more rows than columns).  We use it in Trilinos as an
  /// orthogonalization method for Epetra_MultiVector and
  /// Tpetra::MultiVector.  (In this context, TSQR is provided as an
  /// "OrthoManager" in Anasazi and Belos; you do not have to use it
  /// directly.)  For details, see e.g., our 2008 University of
  /// California Berkeley technical report (Demmel, Grigori, Hoemmen,
  /// and Langou), or our SC09 paper (Demmel, Hoemmen, Mohiyuddin, and
  /// Yelick).
  ///
  /// SequentialTsqr implements the "sequential TSQR" algorithm of the
  /// aforementioned 2008 technical report.  It breaks up the matrix
  /// by rows into "cache blocks," and iterates over consecutive cache
  /// blocks.  The input matrix may be in either the conventional
  /// LAPACK-style column-major layout, or in a "cache-blocked"
  /// layout.  We provide conversion routines between these two
  /// formats.  Users should not attempt to construct a matrix in the
  /// latter format themselves.  In our experience, the performance
  /// difference between the two formats is not significant, but this
  /// may be different on different architectures.
  ///
  /// SequentialTsqr is designed to be used as the "intranode TSQR"
  /// part of the full TSQR implementation in \c Tsqr.  The \c Tsqr
  /// class can use any of various intranode TSQR implementations.
  /// SequentialTsqr is an appropriate choice when running in MPI-only
  /// mode.  Other intranode TSQR implementations, such as \c TbbTsqr,
  /// are appropriate for hybrid parallelism (MPI + threads).
  ///
  /// SequentialTsqr is unlikely to benefit from a multithreaded BLAS
  /// implementation.  In fact, implementations of LAPACK's QR
  /// factorization generally do not show performance benefits from
  /// multithreading when factoring tall skinny matrices.  (See our
  /// SC09 paper and my IPDPS 2011 paper.)  For thread-level
  /// parallelism, use \c TbbTsqr.
  ///
  /// \note To implementers: SequentialTsqr cannot currently be a
  ///   ParameterListAcceptorDefaultBase, because the latter uses RCP,
  ///   and RCPs (more specifically, their reference counts) are not
  ///   currently thread safe.  TbbTsqr uses SequentialTsqr in
  ///   parallel to implement each thread's cache-blocked TSQR.
  ///
  template<class LocalOrdinal, class Scalar>
  class SequentialTsqr : 
    public NodeTsqr<LocalOrdinal, Scalar, std::vector<std::vector<Scalar> > >
  {
  public:
    typedef LocalOrdinal ordinal_type;
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef typename NodeTsqr<LocalOrdinal, Scalar, std::vector<std::vector<Scalar> > >::factor_output_type FactorOutput;

  private:
    typedef typename FactorOutput::const_iterator FactorOutputIter;
    typedef typename FactorOutput::const_reverse_iterator FactorOutputReverseIter;
    typedef MatView<LocalOrdinal, Scalar> mat_view;
    typedef ConstMatView<LocalOrdinal, Scalar> const_mat_view;
    typedef std::pair<mat_view, mat_view> block_pair_type;
    typedef std::pair<const_mat_view, const_mat_view> const_block_pair_type;

    /// \brief Factor the first cache block of the matrix.
    ///
    /// Compute the QR factorization of the first cache block A_top.
    /// Overwrite the upper triangle of A_top with the R factor, and
    /// return a view of the R factor (stored in place in A_top).
    /// Overwrite the (strict) lower triangle of A_top, and the
    /// A_top.ncols() entries of tau, with an implicit representation
    /// of the Q factor.
    ///
    /// \param combine [in/out] Implementation of linear algebra
    ///   primitives.  This is non-const because some Combine
    ///   implementations use scratch space.
    ///
    /// \param A_top [in/out] On input: the first (topmost) cache
    ///   block of the matrix.  Prerequisite: A_top.nrows() >=
    ///   A.top.ncols().  On output, the upper triangle of A_top is
    ///   overwritten with the R factor, and the lower trapezoid of
    ///   A_top is overwritten with part of the implicit
    ///   representation of the Q factor.
    ///
    /// \param tau [out] Array of length >= A_top.ncols().  On output:
    ///   the TAU array (see the LAPACK documentation for _GEQRF).
    ///
    /// \param work [out] Workspace array of length >= A_top.ncols().
    ///
    /// \return A view of the upper triangle of A_top, containing the
    ///   R factor.
    mat_view
    factor_first_block (Combine<LocalOrdinal, Scalar>& combine,
			mat_view& A_top,
			std::vector<Scalar>& tau,
			std::vector<Scalar>& work) const
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
    apply_first_block (Combine<LocalOrdinal, Scalar>& combine,
		       const ApplyType& applyType,
		       const const_mat_view& Q_first,
		       const std::vector<Scalar>& tau,
		       mat_view& C_first,
		       std::vector<Scalar>& work) const
    {
      const LocalOrdinal nrowsLocal = Q_first.nrows();
      combine.apply_first (applyType, nrowsLocal, C_first.ncols(), 
			   Q_first.ncols(), Q_first.get(), Q_first.lda(),
			   &tau[0], C_first.get(), C_first.lda(), &work[0]);
    }

    void
    combine_apply (Combine<LocalOrdinal, Scalar>& combine,
		   const ApplyType& apply_type,
		   const const_mat_view& Q_cur,
		   const std::vector<Scalar>& tau,
		   mat_view& C_top,
		   mat_view& C_cur,
		   std::vector<Scalar>& work) const
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
    combine_factor (Combine<LocalOrdinal, Scalar>& combine,
		    mat_view& R,
		    mat_view& A_cur,
		    std::vector<Scalar>& tau,
		    std::vector<Scalar>& work) const
    {
      const LocalOrdinal nrows_local = A_cur.nrows();
      const LocalOrdinal ncols = A_cur.ncols();

      combine.factor_inner (nrows_local, ncols, R.get(), R.lda(), 
			    A_cur.get(), A_cur.lda(), &tau[0], 
			    &work[0]);
    }

  public:

    /// \brief The standard constructor.
    ///
    /// \param cacheBlockSize [in] Size in bytes of the cache block
    ///   to use in the sequential TSQR factorization.  If 0, the
    ///   implementation will pick a reasonable size.
    ///
    /// \param sizeOfScalar [in] The number of bytes required to store
    ///   a Scalar value.  This is used to compute the dimensions of
    ///   cache blocks.  If sizeof(Scalar) correctly reports the size
    ///   of the representation of Scalar in memory, you can use the
    ///   default.  The default is correct for float, double, and any
    ///   of various fixed-length structs (like double-double and
    ///   quad-double).  It should also work for std::complex<T> where
    ///   T is anything in the previous sentence's list.  It does
    ///   <it>not</it> work for arbitrary-precision types whose
    ///   storage is dynamically allocated, even if the amount of
    ///   storage is a constant.  In the latter case, you should
    ///   specify a nondefault value.  
    ///
    /// \note sizeOfScalar affects performance, not correctness (more
    ///   or less -- it should never be zero, for example).  It's OK
    ///   for it to be a slight overestimate.  Being much too big may
    ///   affect performance by underutilizing the cache.  Being too
    ///   small may also affect performance by thrashing the cache.
    ///
    /// \note If Scalar is an arbitrary-precision type whose
    ///   representation length can change at runtime, you should
    ///   construct a new SequentialTsqr object whenever the
    ///   representation length changes.
    SequentialTsqr (const size_t cacheBlockSize = 0,
		    const size_t sizeOfScalar = sizeof(Scalar)) :
      strategy_ (cacheBlockSize, sizeOfScalar)
    {}

    /// \brief Alternate constructor for a given cache blocking strategy.
    ///
    /// The cache blocking strategy stores the same information as
    /// would be passed into the standard constructor: the cache block
    /// size, and the size of the Scalar type.
    ///
    /// \param strategy [in] Cache blocking strategy to use (copied).
    ///
    SequentialTsqr (const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy) :
      strategy_ (strategy) 
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

    /// \brief Does factor() compute R with nonnegative diagonal?
    ///
    /// See the \c NodeTsqr documentation for details.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      typedef Combine<LocalOrdinal, Scalar> combine_type;
      return combine_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    /// \brief Cache block size in bytes.
    ///
    /// This may be different than the cache block size argument
    /// specified in the constructor.  SequentialTsqr treats that as a
    /// hint, not a command.
    size_t cache_block_size () const { 
      return strategy_.cache_block_size(); 
    }

    /// \brief Compute QR factorization (implicitly stored Q factor) of A.
    ///
    /// Compute the QR factorization in place of the nrows by ncols
    /// matrix A, with nrows >= ncols.  The matrix A is stored either
    /// in column-major order (the default) or with contiguous
    /// column-major cache blocks, with leading dimension lda >=
    /// nrows.  Write the resulting R factor to the top block of A (in
    /// place).  (You can get a view of this via the top_block()
    /// method.)  Everything below the upper triangle of A is
    /// overwritten with part of the implicit representation of the Q
    /// factor.  The other part of that representation is returned.
    ///
    /// \param nrows [in] Number of rows in the matrix A.
    /// \param ncols [in] Number of columns in the matrix A.
    /// \param A [in/out] On input: the nrows by ncols matrix to
    ///   factor.  On output: part of the representation of the
    ///   implicitly stored Q factor.
    /// \param lda [in] Leading dimension of A, if A is stored in
    ///   column-major order.  Otherwise its value doesn't matter.
    /// \param contiguous_cache_blocks [in] Whether the matrix A is
    ///   stored in a contiguously cache-blocked format.
    ///
    /// \return Part of the representation of the implicitly stored Q
    ///   factor.  The complete representation includes A (on output).
    ///   The FactorOutput and A go together.
    FactorOutput
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols,
	    Scalar A[],
	    const LocalOrdinal lda, 
	    const bool contiguous_cache_blocks) const
    {
      CacheBlocker<LocalOrdinal, Scalar> blocker (nrows, ncols, strategy_);
      Combine<LocalOrdinal, Scalar> combine;
      std::vector<Scalar> work (ncols);
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
      std::vector<Scalar> tau_first (ncols);
      mat_view R_view = factor_first_block (combine, A_cur, tau_first, work);
      tau_arrays.push_back (tau_first);

      while (! A_rest.empty())
	{
	  A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	  std::vector<Scalar> tau (ncols);
	  combine_factor (combine, R_view, A_cur, tau, work);
	  tau_arrays.push_back (tau);
	}
      return tau_arrays;
    }

    /// \brief Extract R factor from \c factor() results.
    ///
    /// The five-argument version of \c factor() leaves the R factor
    /// in place in the matrix A.  This method copies the R factor out
    /// of A into a separate matrix R in column-major order
    /// (regardless of whether A was stored with contiguous cache
    /// blocks).
    void
    extract_R (const LocalOrdinal nrows,
	       const LocalOrdinal ncols,
	       const Scalar A[],
	       const LocalOrdinal lda,
	       Scalar R[],
	       const LocalOrdinal ldr,
	       const bool contiguous_cache_blocks) const
    {
      const_mat_view A_view (nrows, ncols, A, lda);

      // Identify top cache block of A
      const_mat_view A_top = top_block (A_view, contiguous_cache_blocks);

      // Fill R (including lower triangle) with zeros.
      fill_matrix (ncols, ncols, R, ldr, Teuchos::ScalarTraits<Scalar>::zero());

      // Copy out the upper triangle of the R factor from A into R.
      copy_upper_triangle (ncols, ncols, R, ldr, A_top.get(), A_top.lda());
    }

    /// \brief Compute the QR factorization of the matrix A.
    ///
    /// See the \c NodeTsqr documentation for details.  This version
    /// of factor() is more useful than the five-argument version,
    /// when using SequentialTsqr as the intranode TSQR implementation
    /// in \c Tsqr.  The five-argument version is more useful when
    /// using SequentialTsqr inside of another intranode TSQR
    /// implementation, such as \c TbbTsqr.
    FactorOutput
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols,
	    Scalar A[],
	    const LocalOrdinal lda, 
	    Scalar R[],
	    const LocalOrdinal ldr,
	    const bool contiguous_cache_blocks) const
    {
      CacheBlocker<LocalOrdinal, Scalar> blocker (nrows, ncols, strategy_);
      Combine<LocalOrdinal, Scalar> combine;
      std::vector<Scalar> work (ncols);
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
      std::vector<Scalar> tau_first (ncols);
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


    /// \brief The number of cache blocks that factor() would use.
    ///
    /// The \c factor() method breaks the input matrix A into one or
    /// more cache blocks.  This method reports how many cache blocks
    /// \c factor() would use, without actually factoring the matrix.
    ///
    /// \param nrows [in] Number of rows in the matrix A.
    /// \param ncols [in] Number of columns in the matrix A.
    /// \param A [in] The matrix A.  If contiguous_cache_blocks is
    ///   false, A is stored in column-major order; otherwise, A is
    ///   stored with contiguous cache blocks (as the \c cache_block()
    ///   method would do).
    /// \param lda [in] If the matrix A is stored in column-major
    ///   order: the leading dimension (a.k.a. stride) of A.
    ///   Otherwise, the value of this parameter doesn't matter.
    /// \param contiguous_cache_blocks [in] Whether the cache blocks
    ///   in the matrix A are stored contiguously.
    ///
    /// \return Number of cache blocks in the matrix A: a positive integer.
    LocalOrdinal
    factor_num_cache_blocks (const LocalOrdinal nrows,
			     const LocalOrdinal ncols,
			     const Scalar A[],
			     const LocalOrdinal lda, 
			     const bool contiguous_cache_blocks) const
    {
      CacheBlocker<LocalOrdinal, Scalar> blocker (nrows, ncols, strategy_);
      LocalOrdinal count = 0;

      const_mat_view A_rest (nrows, ncols, A, lda);
      if (A_rest.empty())
	return count;

      const_mat_view A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
      ++count; // first factor step

      while (! A_rest.empty())
	{
	  A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	  ++count; // next factor step
	}
      return count;
    }

    /// \brief Apply the implicit Q factor to the matrix C.
    ///
    /// See the \c NodeTsqr documentation for details.
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
	   const bool contiguous_cache_blocks) const
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
      CacheBlocker<LocalOrdinal, Scalar> blocker (nrows, ncols_Q, strategy_);
      LAPACK<LocalOrdinal, Scalar> lapack;
      Combine<LocalOrdinal, Scalar> combine;

      const bool transposed = apply_type.transposed();
      const FactorOutput& tau_arrays = factor_output; // rename for encapsulation
      std::vector<Scalar> work (ncols_C);
      
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

      // Identify the top ncols_C by ncols_C block of C.  C_rest is
      // not modified.
      mat_view C_top = blocker.top_block (C_rest, contiguous_cache_blocks);

      if (transposed)
	{
	  const_mat_view Q_cur = blocker.split_top_block (Q_rest, contiguous_cache_blocks);
	  mat_view C_cur = blocker.split_top_block (C_rest, contiguous_cache_blocks);

	  // Apply the topmost block of Q.
	  FactorOutputIter tau_iter = tau_arrays.begin();
	  const std::vector<Scalar>& tau = *tau_iter++;
	  apply_first_block (combine, apply_type, Q_cur, tau, C_cur, work);

	  while (! Q_rest.empty())
	    {
	      Q_cur = blocker.split_top_block (Q_rest, contiguous_cache_blocks);
	      C_cur = blocker.split_top_block (C_rest, contiguous_cache_blocks);
	      combine_apply (combine, apply_type, Q_cur, *tau_iter++, C_top, C_cur, work);
	    }
	}
      else
	{
	  // Start with the last local Q factor and work backwards up the matrix.
	  FactorOutputReverseIter tau_iter = tau_arrays.rbegin();

	  const_mat_view Q_cur = blocker.split_bottom_block (Q_rest, contiguous_cache_blocks);
	  mat_view C_cur = blocker.split_bottom_block (C_rest, contiguous_cache_blocks);

	  while (! Q_rest.empty())
	    {
	      combine_apply (combine, apply_type, Q_cur, *tau_iter++, C_top, C_cur, work);
	      Q_cur = blocker.split_bottom_block (Q_rest, contiguous_cache_blocks);
	      C_cur = blocker.split_bottom_block (C_rest, contiguous_cache_blocks);
	    }
	  // Apply to last (topmost) cache block.
	  apply_first_block (combine, apply_type, Q_cur, *tau_iter++, C_cur, work);
	}
    }

    /// \brief Compute the explicit Q factor from the result of factor().
    ///
    /// See the \c NodeTsqr documentation for details.
    void
    explicit_Q (const LocalOrdinal nrows,
		const LocalOrdinal ncols_Q,
		const Scalar Q[],
		const LocalOrdinal ldq,
		const FactorOutput& factor_output,
		const LocalOrdinal ncols_C,
		Scalar C[],
		const LocalOrdinal ldc,
		const bool contiguous_cache_blocks) const
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

    /// \brief Compute Q := Q*B.
    ///
    /// See the \c NodeTsqr documentation for details.
    void
    Q_times_B (const LocalOrdinal nrows,
	       const LocalOrdinal ncols,
	       Scalar Q[],
	       const LocalOrdinal ldq,
	       const Scalar B[],
	       const LocalOrdinal ldb,
	       const bool contiguous_cache_blocks) const
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
      CacheBlocker<LocalOrdinal, Scalar> blocker (nrows, ncols, strategy_);
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
      CacheBlocker<LocalOrdinal, Scalar> blocker (nrows, ncols, strategy_);
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
		     const bool contiguous_cache_blocks) const
    {
      CacheBlocker<LocalOrdinal, Scalar> blocker (nrows, ncols, strategy_);
      blocker.fill_with_zeros (nrows, ncols, A, lda, contiguous_cache_blocks);
    }

  protected:

    /// \brief Return the topmost cache block of the matrix C.
    ///
    /// NodeTsqr's top_block() method must be implemented using
    /// subclasses' const_top_block() method, since top_block() is a
    /// template method and template methods cannot be virtual.
    ///
    /// \param C [in] View of a matrix, with at least as many rows as
    ///   columns.
    /// \param contiguous_cache_blocks [in] Whether the cache blocks
    ///   of C are stored contiguously.
    ///
    /// \return View of the topmost cache block of the matrix C.
    ConstMatView<LocalOrdinal, Scalar>
    const_top_block (const ConstMatView<LocalOrdinal, Scalar>& C, 
		     const bool contiguous_cache_blocks) const 
    {
      // The CacheBlocker object knows how to construct a view of the
      // top cache block of C.  This is complicated because cache
      // blocks (in C) may or may not be stored contiguously.  If they
      // are stored contiguously, the CacheBlocker knows the right
      // layout, based on the cache blocking strategy.
      typedef CacheBlocker<LocalOrdinal, Scalar> blocker_type;
      blocker_type blocker (C.nrows(), C.ncols(), strategy_);

      // C_top_block is a view of the topmost cache block of C.
      // C_top_block should have >= ncols rows, otherwise either cache
      // blocking is broken or the input matrix C itself had fewer
      // rows than columns.
      const_mat_view C_top_block = 
	blocker.top_block (C, contiguous_cache_blocks);
      return C_top_block;
    }

  private:
    //! Strategy object that helps us cache block matrices.
    CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
  };
  
} // namespace TSQR

#endif // __TSQR_Tsqr_SequentialTsqr_hpp
