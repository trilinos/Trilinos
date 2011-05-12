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

#ifndef __TSQR_TbbTsqr_hpp
#define __TSQR_TbbTsqr_hpp

#include <TbbTsqr_TbbParallelTsqr.hpp>
#include <Tsqr_TimeStats.hpp>
#include <Teuchos_Time.hpp>
// #include <TbbRecursiveTsqr.hpp>

#include <stdexcept>
#include <string>
#include <utility> // std::pair
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {

    /// \class TbbTsqr
    /// \brief Intranode TSQR, parallelized with Intel TBB
    ///
    /// TSQR factorization for a dense, tall and skinny matrix stored
    /// on a single node.  Parallelized using Intel's Threading
    /// Building Blocks.
    ///
    /// \note TSQR only needs to know about the local ordinal type
    ///   (LocalOrdinal), not about the global ordinal type.
    ///   TimerType may be any class with the same interface as
    ///   TrivialTimer; it times the divide-and-conquer base cases
    ///   (the operations on each CPU core within the thread-parallel
    ///   implementation).
    template< class LocalOrdinal, class Scalar, class TimerType = Teuchos::Time >
    class TbbTsqr : public Teuchos::Describable {
    private:
      /// \brief Implementation of TBB TSQR.
      ///
      /// If you don't have TBB available, you can test this class by
      /// substituting in a TbbRecursiveTsqr<LocalOrdinal, Scalar>
      /// object.  That is a nonparallel implementation that emulates
      /// the control flow of TbbParallelTsqr.  If you do this, you
      /// should also change the FactorOutput public typedef.
      ///
      /// \note This is NOT a use of the pImpl idiom.  
      TbbParallelTsqr< LocalOrdinal, Scalar, TimerType > impl_;

      // Collected running statistcs on various computations
      mutable TimeStats factorStats_, applyStats_, explicitQStats_, cacheBlockStats_, unCacheBlockStats_;

      // Timers for various computations
      mutable TimerType factorTimer_, applyTimer_, explicitQTimer_, cacheBlockTimer_, unCacheBlockTimer_;

    public:
      typedef Scalar scalar_type;
      typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
      typedef LocalOrdinal ordinal_type;

      /// \typedef FactorOutput
      /// \brief Type of partial output of TBB TSQR.
      ///
      /// If you don't have TBB available, you can test this class by
      /// substituting in "typename TbbRecursiveTsqr<LocalOrdinal,
      /// Scalar>::FactorOutput" for the typedef's definition.  If you
      /// do this, you should also change the type of \c impl_ above.
      typedef typename TbbParallelTsqr<LocalOrdinal, Scalar, TimerType>::FactorOutput FactorOutput;

      //! (Max) number of cores used for the factorization.
      size_t ncores() const { return impl_.ncores(); }

      //! Cache size hint (in bytes) used for the factorization.
      size_t cache_size_hint() const { return impl_.cache_size_hint(); }

      /// \brief Cache size hint (in bytes) used for the factorization.
      ///
      /// This method is deprecated, because the name is misleading.
      /// Please call \c cache_size_hint() instead.
      size_t cache_block_size() const { return impl_.cache_size_hint(); }

      /// \brief Constructor; sets up tuning parameters.
      ///
      /// \param numCores [in] Maximum number of processing cores to use
      ///   when factoring the matrix.  Fewer cores may be used if the
      ///   matrix is not big enough to justify their use.
      ///
      /// \param cacheSizeHint [in] Cache block size hint (in bytes)
      ///   to use in the sequential part of TSQR.  If zero or not
      ///   specified, a reasonable default is used.  If each CPU core
      ///   has a private cache, that cache's size (minus a little
      ///   wiggle room) would be the appropriate value for this
      ///   parameter.  Set to zero for the implementation to choose a
      ///   reasonable default.
      TbbTsqr (const size_t numCores,
	       const size_t cacheSizeHint = 0) :
	impl_ (numCores, cacheSizeHint),
	factorTimer_ ("TbbTsqr::factor"),
	applyTimer_ ("TbbTsqr::apply"),
	explicitQTimer_ ("TbbTsqr::explicit_Q"),
	cacheBlockTimer_ ("TbbTsqr::cache_block"),
	unCacheBlockTimer_ ("TbbTsqr::un_cache_block")
      {}

      /// Whether or not this QR factorization produces an R factor
      /// with all nonnegative diagonal entries.
      static bool QR_produces_R_factor_with_nonnegative_diagonal() {
	typedef TbbParallelTsqr< LocalOrdinal, Scalar, TimerType > impl_type;
	return impl_type::QR_produces_R_factor_with_nonnegative_diagonal();
      }

      /// \brief One-line description of this object.
      ///
      /// This implements Teuchos::Describable::description().  For now,
      /// SequentialTsqr uses the default implementation of
      /// Teuchos::Describable::describe().
      std::string description () const {
	using std::endl;

	// SequentialTsqr also implements Describable, so if you
	// decide to implement describe(), you could call
	// SequentialTsqr's describe() and get a nice hierarchy of
	// descriptions.
	std::ostringstream os;
	os << "Intranode Tall Skinny QR (TSQR): "
	   << "Intel Threading Building Blocks (TBB) implementation"
	   << ", max " << ncores() << "-way parallelism"
	   << ", cache size hint of " << cache_size_hint() << " bytes.";
	return os.str();
      }

      void
      cache_block (const LocalOrdinal nrows,
		   const LocalOrdinal ncols, 
		   Scalar A_out[],
		   const Scalar A_in[],
		   const LocalOrdinal lda_in) const
      {
	cacheBlockTimer_.start(true);
	impl_.cache_block (nrows, ncols, A_out, A_in, lda_in);
	cacheBlockStats_.update (cacheBlockTimer_.stop());
      }

      void
      un_cache_block (const LocalOrdinal nrows,
		      const LocalOrdinal ncols,
		      Scalar A_out[],
		      const LocalOrdinal lda_out,		    
		      const Scalar A_in[]) const
      {
	unCacheBlockTimer_.start(true);
	impl_.un_cache_block (nrows, ncols, A_out, lda_out, A_in);
	unCacheBlockStats_.update (unCacheBlockTimer_.stop());
      }

      void
      fill_with_zeros (const LocalOrdinal nrows,
		       const LocalOrdinal ncols,
		       Scalar C[],
		       const LocalOrdinal ldc, 
		       const bool contiguous_cache_blocks) const
      {
	impl_.fill_with_zeros (nrows, ncols, C, ldc, contiguous_cache_blocks);
      }

      template< class MatrixViewType >
      MatrixViewType
      top_block (const MatrixViewType& C, 
		 const bool contiguous_cache_blocks) const
      {
	return impl_.top_block (C, contiguous_cache_blocks);
      }

      /// \brief Compute QR factorization of the dense matrix A
      ///
      /// Compute the QR factorization of the dense matrix A.
      ///
      /// \param nrows [in] Number of rows of A.  
      ///   Precondition: nrows >= ncols.
      ///
      /// \param ncols [in] Number of columns of A.
      ///   Precondition: nrows >= ncols.
      ///
      /// \param A [in,out] On input, the matrix to factor, stored as a
      ///   general dense matrix in column-major order.  On output,
      ///   overwritten with an implicit representation of the Q factor.
      ///
      /// \param lda [in] Leading dimension of A.  
      ///   Precondition: lda >= nrows.
      ///
      /// \param R [out] The final R factor of the QR factorization of
      ///   the matrix A.  An ncols by ncols upper triangular matrix
      ///   stored in column-major order, with leading dimension ldr.
      ///
      /// \param ldr [in] Leading dimension of the matrix R.
      ///
      /// \param b_contiguous_cache_blocks [in] Whether cache blocks are
      ///   stored contiguously in the input matrix A and the output
      ///   matrix Q (of explicit_Q()).  If not and you want them to be,
      ///   you should use the cache_block() method to copy them into
      ///   that format.  You may use the un_cache_block() method to
      ///   copy them out of that format into the usual column-oriented
      ///   format.
      ///
      /// \return FactorOutput struct, which together with the data in A
      ///   form an implicit representation of the Q factor.  They
      ///   should be passed into the apply() and explicit_Q() functions
      ///   as the "factor_output" parameter.
      FactorOutput
      factor (const LocalOrdinal nrows,
	      const LocalOrdinal ncols, 
	      Scalar A[],
	      const LocalOrdinal lda,
	      Scalar R[],
	      const LocalOrdinal ldr,
	      const bool contiguous_cache_blocks) const
      {
	factorTimer_.start(true);
	return impl_.factor (nrows, ncols, A, lda, R, ldr, contiguous_cache_blocks);
	factorStats_.update (factorTimer_.stop());
      }

      /// \brief Apply Q factor to the global dense matrix C
      ///
      /// Apply the Q factor (computed by factor() and represented
      /// implicitly) to the dense matrix C.
      ///
      /// \param apply_type [in] Whether to compute Q*C, Q^T * C, or
      ///   Q^H * C.
      ///
      /// \param nrows [in] Number of rows of the matrix C and the
      ///   matrix Q.  Precondition: nrows >= ncols_Q, ncols_C.
      ///
      /// \param ncols_Q [in] Number of columns of Q
      ///
      /// \param Q [in] Same as the "A" output of factor()
      ///
      /// \param ldq [in] Same as the "lda" input of factor()
      ///
      /// \param factor_output [in] Return value of factor()
      ///
      /// \param ncols_C [in] Number of columns in C.
      ///   Precondition: nrows_local >= ncols_C.
      ///
      /// \param C [in,out] On input, the matrix C, stored as a general
      ///   dense matrix in column-major order.  On output, overwritten
      ///   with op(Q)*C, where op(Q) = Q or Q^T.
      ///
      /// \param ldc [in] Leading dimension of C.  
      ///   Precondition: ldc_local >= nrows_local.
      ///   Not applicable if C is cache-blocked in place.
      ///
      /// \param contiguous_cache_blocks [in] Whether or not cache
      ///   blocks of Q and C are stored contiguously (default:
      ///   false).
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
	applyTimer_.start(true);
	impl_.apply (apply_type, nrows, ncols_Q, Q, ldq, factor_output, 
		     ncols_C, C, ldc, contiguous_cache_blocks);
	applyStats_.update (applyTimer_.stop());
      }

      /// \brief Compute the explicit Q factor from factor()
      ///
      /// Compute the explicit version of the Q factor computed by
      /// factor() and represented implicitly (via Q_in and
      /// factor_output).
      ///
      /// \param nrows [in] Number of rows of the matrix Q_in.  Also,
      ///   the number of rows of the output matrix Q_out.
      ///   Precondition: nrows >= ncols_Q_in.
      ///
      /// \param ncols_Q_in [in] Number of columns in the original matrix
      ///   A, whose explicit Q factor we are computing.
      ///   Precondition: nrows >= ncols_Q_in.
      ///
      /// \param Q_local_in [in] Same as A output of factor().
      ///
      /// \param ldq_local_in [in] Same as lda input of factor()
      ///
      /// \param ncols_Q_out [in] Number of columns of the explicit Q
      ///   factor to compute.
      ///
      /// \param Q_out [out] The explicit representation of the Q factor.
      ///
      /// \param ldq_out [in] Leading dimension of Q_out.
      ///
      /// \param factor_output [in] Return value of factor().
      void
      explicit_Q (const LocalOrdinal nrows,
		  const LocalOrdinal ncols_Q_in,
		  const Scalar Q_in[],
		  const LocalOrdinal ldq_in,
		  const FactorOutput& factor_output,
		  const LocalOrdinal ncols_Q_out,
		  Scalar Q_out[],
		  const LocalOrdinal ldq_out,
		  const bool contiguous_cache_blocks) const
      {
	explicitQTimer_.start(true);
	impl_.explicit_Q (nrows, ncols_Q_in, Q_in, ldq_in, factor_output,
			  ncols_Q_out, Q_out, ldq_out, contiguous_cache_blocks);
	explicitQStats_.update (explicitQTimer_.stop());
      }

      /// \brief Compute Q*B
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
		 const bool contiguous_cache_blocks) const
      {
	impl_.Q_times_B (nrows, ncols, Q, ldq, B, ldb, contiguous_cache_blocks);
      }

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
	return impl_.reveal_R_rank (ncols, R, ldr, U, ldu, tol);
      }

      /// \brief Rank-revealing decomposition
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
		   const bool contiguous_cache_blocks) const
      {
	return impl_.reveal_rank (nrows, ncols, Q, ldq, R, ldr, tol, 
				  contiguous_cache_blocks);
      }

      double
      min_seq_factor_timing () const { return impl_.min_seq_factor_timing(); }
      double
      max_seq_factor_timing () const { return impl_.max_seq_factor_timing(); }
      double
      min_seq_apply_timing () const { return impl_.min_seq_apply_timing(); }
      double
      max_seq_apply_timing () const { return impl_.max_seq_apply_timing(); }

      void getStats (std::vector< TimeStats >& stats) {
	const int numStats = 5;
	stats.resize (numStats);
	stats[0] = factorStats_;
	stats[1] = applyStats_;
	stats[2] = explicitQStats_;
	stats[3] = cacheBlockStats_;
	stats[4] = unCacheBlockStats_;
      }

      void getStatsLabels (std::vector< std::string >& labels) {
	const int numStats = 5;
	labels.resize (numStats);
	labels[0] = factorTimer_.name();
	labels[1] = applyTimer_.name();
	labels[2] = explicitQTimer_.name();
	labels[3] = cacheBlockTimer_.name();
	labels[4] = unCacheBlockTimer_.name();
      }

    }; // class TbbTsqr

  } // namespace TBB
} // namespace TSQR

#endif // __TSQR_TbbTsqr_hpp
