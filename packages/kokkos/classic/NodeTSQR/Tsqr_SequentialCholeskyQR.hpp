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

#ifndef __TSQR_Tsqr_SequentialCholeskyQR_hpp
#define __TSQR_Tsqr_SequentialCholeskyQR_hpp

#include <Tsqr_MatView.hpp>
#include <Tsqr_Blas.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_CacheBlockingStrategy.hpp>
#include <Tsqr_CacheBlocker.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_Util.hpp>

#include <string>
#include <utility>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class SequentialCholeskyQR
  /// \brief Cache-blocked sequential implementation of CholeskyQR.
  ///
  /// CholeskyQR works like this: given an input matrix A with no
  /// fewer rows than columns,
  /// - Compute the Gram matrix of A: \f$H = A^* A\f$
  /// - Compute the (upper triangular) Cholesky factorization of H:
  ///   \f$H = R^* R\f$
  /// - Compute \f$Q = A R^{-1}\f$
  template<class LocalOrdinal, class Scalar>
  class SequentialCholeskyQR {
  private:
    typedef MatView< LocalOrdinal, Scalar > mat_view;
    typedef ConstMatView< LocalOrdinal, Scalar > const_mat_view;

  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;

    /// \typedef FactorOutput
    /// \brief Return value of \c factor().
    ///
    /// Here, FactorOutput is just a minimal object whose value is
    /// irrelevant, so that this class' interface looks like that of
    /// \c SequentialTsqr.
    typedef int FactorOutput;

    /// \brief Cache size hint (in bytes).
    ///
    /// This method is deprecated, because the name is misleading.
    /// Please call \c cache_size_hint() instead.
    size_t TEUCHOS_DEPRECATED cache_block_size () const { 
      return strategy_.cache_size_hint(); 
    }

    //! Cache size hint (in bytes).
    size_t cache_size_hint () const { return strategy_.cache_size_hint(); }

    /// \brief Constructor
    ///
    /// \param theCacheSizeHint [in] Cache size hint in bytes.  If 0,
    ///   the implementation will pick a reasonable size, which may be
    ///   queried by calling cache_size_hint().
    SequentialCholeskyQR (const size_t theCacheSizeHint = 0) :
      strategy_ (theCacheSizeHint)
    {}

    /// \brief Whether the R factor has a nonnegative diagonal.
    /// 
    /// The \c factor() method computes a QR factorization of the
    /// input matrix A.  Some, but not all methods for computing a QR
    /// factorization produce an R factor with a nonnegative diagonal.
    /// This class' implementation does, because the R factor comes
    /// from a Cholesky factorization.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return true;
    }

    /// \brief Compute the QR factorization of the matrix A.
    ///
    /// Compute the QR factorization of the nrows by ncols matrix A,
    /// with nrows >= ncols, stored either in column-major order (the
    /// default) or as contiguous column-major cache blocks, with
    /// leading dimension lda >= nrows.
    FactorOutput
    factor (const LocalOrdinal nrows,
	    const LocalOrdinal ncols,
	    const Scalar A[],
	    const LocalOrdinal lda, 
	    Scalar R[],
	    const LocalOrdinal ldr,
	    const bool contiguous_cache_blocks = false)
    {
      CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
      LAPACK< LocalOrdinal, Scalar > lapack;
      BLAS< LocalOrdinal, Scalar > blas;
      std::vector< Scalar > work (ncols);
      Matrix< LocalOrdinal, Scalar > ATA (ncols, ncols, Scalar(0));
      FactorOutput retval (0);

      if (contiguous_cache_blocks)
	{
	  // Compute ATA := A^T * A, by iterating through the cache
	  // blocks of A from top to bottom.
	  //
	  // We say "A_rest" because it points to the remaining part of
	  // the matrix left to process; at the beginning, the "remaining"
	  // part is the whole matrix, but that will change as the
	  // algorithm progresses.
	  mat_view A_rest (nrows, ncols, A, lda);
	  // This call modifies A_rest (but not the actual matrix
	  // entries; just the dimensions and current position).
	  mat_view A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	  // Process the first cache block: ATA := A_cur^T * A_cur
	  blas.GEMM ("T", "N", ncols, ncols, A_cur.nrows(), 
		     Scalar(1), A_cur.get(), A_cur.lda(), A_cur.get(), A_cur.lda(),
		     Scalar(0), ATA.get(), ATA.lda());
	  // Process the remaining cache blocks in order.
	  while (! A_rest.empty())
	    {
	      A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	      // ATA := ATA + A_cur^T * A_cur
	      blas.GEMM ("T", "N", ncols, ncols, A_cur.nrows(), 
			 Scalar(1), A_cur.get(), A_cur.lda(), A_cur.get(), A_cur.lda(),
			 Scalar(1), ATA.get(), ATA.lda());
	    }
	}
      else
	// Compute ATA := A^T * A, using a single BLAS call.
	blas.GEMM ("T", "N", ncols, ncols, nrows, 
		   Scalar(1), A, lda, A, lda,
		   Scalar(0), ATA.get(), ATA.lda());

      // Compute the Cholesky factorization of ATA in place, so that
      // A^T * A = R^T * R, where R is ncols by ncols upper
      // triangular.
      int info = 0;
      lapack.POTRF ("U", ncols, ATA.get(), ATA.lda(), &info);
      // FIXME (mfh 22 June 2010) The right thing to do here would be
      // to resort to a rank-revealing factorization, as Stathopoulos
      // and Wu (2002) do with their CholeskyQR + symmetric
      // eigensolver factorization.
      if (info != 0)
	throw std::runtime_error("Cholesky factorization failed");

      // Copy out the R factor
      fill_matrix (ncols, ncols, R, ldr, Scalar(0));
      copy_upper_triangle (ncols, ncols, R, ldr, ATA.get(), ATA.lda());

      // Compute A := A * R^{-1}.  We do this in place in A, using
      // BLAS' TRSM with the R factor (form POTRF) stored in the upper
      // triangle of ATA.
      {
	mat_view A_rest (nrows, ncols, A, lda);
	// This call modifies A_rest.
	mat_view A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);

	// Compute A_cur / R (Matlab notation for A_cur * R^{-1}) in place.
	blas.TRSM ("R", "U", "N", "N", A_cur.nrows(), ncols, 
		   Scalar(1), ATA.get(), ATA.lda(), A_cur.get(), A_cur.lda());

	// Process the remaining cache blocks in order.
	while (! A_rest.empty())
	  {
	    A_cur = blocker.split_top_block (A_rest, contiguous_cache_blocks);
	    blas.TRSM ("R", "U", "N", "N", A_cur.nrows(), ncols, 
		       Scalar(1), ATA.get(), ATA.lda(), A_cur.get(), A_cur.lda());
	  }
      }

      return retval;
    }

    /// \param factor_output [in] Not used; just here to match the
    ///   interface of SequentialTsqr.
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
      if (ncols_Q != ncols_C)
	throw std::logic_error("SequentialCholeskyQR::explicit_Q() "
			       "does not work if ncols_C != ncols_Q");
      const LocalOrdinal ncols = ncols_Q;

      if (contiguous_cache_blocks)
	{
	  CacheBlocker< LocalOrdinal, Scalar > blocker (nrows, ncols, strategy_);
	  mat_view C_rest (nrows, ncols, C, ldc);
	  const_mat_view Q_rest (nrows, ncols, Q, ldq);

	  mat_view C_cur = blocker.split_top_block (C_rest, contiguous_cache_blocks);
	  const_mat_view Q_cur = blocker.split_top_block (Q_rest, contiguous_cache_blocks);

	  while (! C_rest.empty())
	    Q_cur.copy (C_cur);
	}
      else
	{
	  mat_view C_view (nrows, ncols, C, ldc);
	  C_view.copy (const_mat_view (nrows, ncols, Q, ldq));
	}
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

#endif // __TSQR_Tsqr_SequentialCholeskyQR_hpp
