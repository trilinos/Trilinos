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

#ifndef __TSQR_TbbTsqr_hpp
#define __TSQR_TbbTsqr_hpp

#include <TbbTsqr_TbbParallelTsqr.hpp>
// #include <TbbRecursiveTsqr.hpp>

#include <stdexcept>
#include <string>
#include <utility> // std::pair
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {

    /// \class TrivialTimer
    /// \brief Satisfies TimerType concept trivially
    class TrivialTimer {
    public:
      TrivialTimer (const std::string& name = std::string("NoName")) {}
      void start() {}
      double stop() { return double(0); }
    };

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
    template< class LocalOrdinal, class Scalar, class TimerType = TrivialTimer >
    class TbbTsqr {
    private:
      // Note: this is NOT a use of the pImpl idiom.  TbbRecursiveTsqr
      // is a nonparallel implementation that emulates the control
      // flow of the parallel implementation TbbParallelTsqr.  The
      // latter depends on the Intel Threading Building Blocks
      // library.
      //
      //TbbRecursiveTsqr< LocalOrdinal, Scalar > impl_;
      TbbParallelTsqr< LocalOrdinal, Scalar, TimerType > impl_;

    public:
      typedef Scalar scalar_type;
      typedef LocalOrdinal ordinal_type;
      // typedef typename TbbRecursiveTsqr< LocalOrdinal, Scalar >::FactorOutput FactorOutput;
      typedef typename TbbParallelTsqr< LocalOrdinal, Scalar, TimerType >::FactorOutput FactorOutput;

      /// (Max) number of cores used for the factorization.
      size_t ncores() const { return impl_.ncores(); }

      /// Cache block size (in bytes) used for the factorization.
      size_t cache_block_size() const { return impl_.cache_block_size(); }

      /// Constructor; sets up tuning parameters
      ///
      /// \param numCores [in] Maximum number of processing cores to use
      ///   when factoring the matrix.  Fewer cores may be used if the
      ///   matrix is not big enough to justify their use.
      /// \param cacheBlockSize [in] Size (in bytes) of cache block to
      ///   use in the sequential part of TSQR.  If zero or not
      ///   specified, a reasonable default is used.  If each core has a
      ///   private cache, that cache's size (minus a little wiggle
      ///   room) would be the appropriate value for this parameter.
      ///   Set to zero for the implementation to choose a default,
      ///   which may or may not give good performance on your platform.
      TbbTsqr (const size_t numCores,
	       const size_t cacheBlockSize = 0) :
	impl_ (numCores, cacheBlockSize)
      {}

      void
      cache_block (const LocalOrdinal nrows,
		   const LocalOrdinal ncols, 
		   Scalar A_out[],
		   const Scalar A_in[],
		   const LocalOrdinal lda_in) const
      {
	impl_.cache_block (nrows, ncols, A_out, A_in, lda_in);
      }

      void
      un_cache_block (const LocalOrdinal nrows,
		      const LocalOrdinal ncols,
		      Scalar A_out[],
		      const LocalOrdinal lda_out,		    
		      const Scalar A_in[]) const
      {
	impl_.un_cache_block (nrows, ncols, A_out, lda_out, A_in);
      }

      void
      fill_with_zeros (const LocalOrdinal nrows,
		       const LocalOrdinal ncols,
		       Scalar C[],
		       const LocalOrdinal ldc, 
		       const bool contiguous_cache_blocks = false) const
      {
	impl_.fill_with_zeros (nrows, ncols, C, ldc, contiguous_cache_blocks);
      }

      template< class MatrixViewType >
      MatrixViewType
      top_block (const MatrixViewType& C, 
		 const bool contiguous_cache_blocks = false) const
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
	      const bool contiguous_cache_blocks = false)
      {
	return impl_.factor (nrows, ncols, A, lda, R, ldr, contiguous_cache_blocks);
      }

      /// \brief Apply Q factor to the global dense matrix C
      ///
      /// Apply the Q factor (computed by factor() and represented
      /// implicitly) to the dense matrix C.
      ///
      /// \param [in] If "N", compute Q*C.  If "T", compute Q^T * C.
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
      apply (const std::string& op,
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
	impl_.apply (op, nrows, ncols_Q, Q, ldq, factor_output, 
		     ncols_C, C, ldc, contiguous_cache_blocks);
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
		  const bool contiguous_cache_blocks = false)
      {
	impl_.explicit_Q (nrows, ncols_Q_in, Q_in, ldq_in, factor_output,
			  ncols_Q_out, Q_out, ldq_out, contiguous_cache_blocks);
      }

      double
      min_seq_factor_timing () const { return impl_.min_seq_factor_timing(); }
      double
      max_seq_factor_timing () const { return impl_.max_seq_factor_timing(); }
      double
      min_seq_apply_timing () const { return impl_.min_seq_apply_timing(); }
      double
      max_seq_apply_timing () const { return impl_.max_seq_apply_timing(); }

    }; // class TbbTsqr

  } // namespace TBB
} // namespace TSQR

#endif // __TSQR_TbbTsqr_hpp
