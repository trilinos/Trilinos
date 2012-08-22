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

#ifndef __TSQR_Tsqr_LocalVerify_hpp
#define __TSQR_Tsqr_LocalVerify_hpp

#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_Blas.hpp>
#include <Tsqr_Util.hpp>

#include <cmath>
#include <limits>
#include <utility> // std::pair, std::make_pair
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template< class Ordinal, class Scalar >
  typename ScalarTraits< Scalar >::magnitude_type
  local_frobenius_norm (const Ordinal nrows_local, 
			const Ordinal ncols,
			const Scalar  A_local[],
			const Ordinal lda_local)
  {
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
    
    // FIXME (mfh 22 Apr 2010) This function does no scaling of
    // intermediate quantities, so it might overflow unnecessarily.
    magnitude_type result (0);
    for (Ordinal j = 0; j < ncols; j++)
      {
	const Scalar* const cur_col = &A_local[j*lda_local];
	for (Ordinal i = 0; i < nrows_local; ++i)
	  {
	    const magnitude_type abs_xi = ScalarTraits< Scalar >::abs (cur_col[i]);
	    result = result + abs_xi * abs_xi;
	  }
      }
    return sqrt (result);
  }


  template< class Ordinal, class Scalar >
  bool
  NaN_in_matrix (const Ordinal nrows, 
		 const Ordinal ncols, 
		 const Scalar A[], 
		 const Ordinal lda)
  {
    // Testing whether a NaN is present in A only makes sense if it is
    // possible for NaNs not to signal.  Otherwise the NaNs would have
    // signalled and we wouldn't need to be here.  Of course perhaps
    // one could change the signal state at runtime, but has_quiet_NaN
    // refers to the possibility of quiet NaNs being able to exist at
    // all.
    if (std::numeric_limits<Scalar>::has_quiet_NaN)
      {
	for (Ordinal j = 0; j < ncols; j++)
	  for (Ordinal i = 0; i < nrows; i++)
	    {
	      if (std::isnan (A[i + j*lda]))
		return true;
	    }
	return false;
      }
    else
      return false;
  }


  template< class Ordinal, class Scalar >
  bool
  NaN_in_matrix (const Ordinal nrows, 
		 const Ordinal ncols, 
		 const std::vector<Scalar>& A, 
		 const Ordinal lda)
  {
    const Scalar* const A_ptr = &A[0];
    return NaN_in_matrix (nrows, ncols, A_ptr, lda);
  }



  template< class Ordinal, class Scalar >
  typename ScalarTraits< Scalar >::magnitude_type
  localOrthogonality (const Ordinal nrows,
		      const Ordinal ncols,
		      const Scalar Q[],
		      const Ordinal ldq)
  {
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
    const Scalar ZERO (0);
    const Scalar ONE (1);

    BLAS<Ordinal, Scalar> blas;
  
    std::vector<Scalar> AbsOrthog (ncols * ncols, std::numeric_limits<Scalar>::quiet_NaN());
    const Ordinal AbsOrthog_stride = ncols;

    // Compute AbsOrthog := Q' * Q - I.  First, compute Q' * Q:
    if (ScalarTraits< Scalar >::is_complex)
      blas.GEMM ("C", "N", ncols, ncols, nrows, ONE, Q, ldq, Q, ldq, 
		 ZERO, &AbsOrthog[0], AbsOrthog_stride);
    else
      blas.GEMM ("T", "N", ncols, ncols, nrows, ONE, Q, ldq, Q, ldq, 
		 ZERO, &AbsOrthog[0], AbsOrthog_stride);

    // Now, compute (Q^T*Q) - I.
    for (Ordinal j = 0; j < ncols; j++)
      AbsOrthog[j + j*AbsOrthog_stride] = AbsOrthog[j + j*AbsOrthog_stride] - ONE;

    // Now AbsOrthog == Q^T * Q - I.  Compute and return its Frobenius norm.
    return local_frobenius_norm (ncols, ncols, &AbsOrthog[0], AbsOrthog_stride);
  }


  
  template< class Ordinal, class Scalar >
  typename ScalarTraits< Scalar >::magnitude_type
  local_relative_orthogonality (const Ordinal nrows,
				const Ordinal ncols,
				const Scalar Q[],
				const Ordinal ldq,
				const typename ScalarTraits< Scalar >::magnitude_type A_norm_F)
  {
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
    const Scalar ZERO (0);
    const Scalar ONE (1);

    const bool relative = false; // whether to scale $\|I-Q^T*Q\|_F$ by $\|A\|_F$
    BLAS<Ordinal, Scalar> blas;
  
    std::vector<Scalar> AbsOrthog (ncols * ncols, std::numeric_limits<Scalar>::quiet_NaN());
    const Ordinal AbsOrthog_stride = ncols;

    // Compute AbsOrthog := Q' * Q - I.  First, compute Q' * Q:
    if (ScalarTraits< Scalar >::is_complex)
      blas.GEMM ("C", "N", ncols, ncols, nrows, ONE, Q, ldq, Q, ldq, 
		 ZERO, &AbsOrthog[0], AbsOrthog_stride);
    else
      blas.GEMM ("T", "N", ncols, ncols, nrows, ONE, Q, ldq, Q, ldq, 
		 ZERO, &AbsOrthog[0], AbsOrthog_stride);

    // Now, compute (Q^T*Q) - I.
    for (Ordinal j = 0; j < ncols; j++)
      AbsOrthog[j + j*AbsOrthog_stride] = AbsOrthog[j + j*AbsOrthog_stride] - ONE;

    // Now AbsOrthog == Q^T * Q - I.  Compute its Frobenius norm.
    const magnitude_type AbsOrthog_norm_F = 
      local_frobenius_norm (ncols, ncols, &AbsOrthog[0], AbsOrthog_stride);

    // Return the orthogonality measure
    return relative ? (AbsOrthog_norm_F / A_norm_F) : AbsOrthog_norm_F;
  }


  template< class Ordinal, class Scalar >
  typename ScalarTraits< Scalar >::magnitude_type
  localResidual (const Ordinal nrows, 
		 const Ordinal ncols,
		 const Scalar A[],
		 const Ordinal lda,
		 const Scalar Q[],
		 const Ordinal ldq,
		 const Scalar R[],
		 const Ordinal ldr)
  {
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

    std::vector< Scalar > AbsResid (nrows * ncols, std::numeric_limits< Scalar >::quiet_NaN());
    const Ordinal AbsResid_stride = nrows;
    BLAS< Ordinal, Scalar > blas;
    const magnitude_type ONE (1);

    // A_copy := A_copy - Q * R
    copy_matrix (nrows, ncols, &AbsResid[0], AbsResid_stride, A, lda);
    blas.GEMM ("N", "N", nrows, ncols, ncols, -ONE, Q, ldq, R, ldr, 
	       ONE, &AbsResid[0], AbsResid_stride);

    return local_frobenius_norm (nrows, ncols, &AbsResid[0], AbsResid_stride);
  }


  template< class Ordinal, class Scalar >
  typename ScalarTraits< Scalar >::magnitude_type
  local_relative_residual (const Ordinal nrows, 
			   const Ordinal ncols,
			   const Scalar A[],
			   const Ordinal lda,
			   const Scalar Q[],
			   const Ordinal ldq,
			   const Scalar R[],
			   const Ordinal ldr,
			   const typename ScalarTraits< Scalar >::magnitude_type A_norm_F)
  {
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

    std::vector< Scalar > AbsResid (nrows * ncols, std::numeric_limits< Scalar >::quiet_NaN());
    const Ordinal AbsResid_stride = nrows;
    BLAS< Ordinal, Scalar > blas;
    const magnitude_type ONE (1);

    // if (b_debug)
    //   cerr << "relative_residual:" << endl;
    // if (matrix_contains_nan (nrows, ncols, A, lda))
    //   cerr << "relative_residual: matrix A contains a NaN" << endl;
    // if (matrix_contains_nan (nrows, ncols, Q, ldq))
    //   cerr << "relative_residual: matrix Q contains a NaN" << endl;
    // if (matrix_contains_nan (ncols, ncols, R, ldr))
    //   cerr << "relative_residual: matrix R contains a NaN" << endl;

    // A_copy := A_copy - Q * R
    copy_matrix (nrows, ncols, &AbsResid[0], AbsResid_stride, A, lda);

    // if (NaN_in_matrix (nrows, ncols, AbsResid, AbsResid_stride))
    //   cerr << "relative_residual: matrix AbsResid := A contains a NaN" << endl;

    blas.GEMM ("N", "N", nrows, ncols, ncols, -ONE, Q, ldq, R, ldr, 
	       ONE, &AbsResid[0], AbsResid_stride);

    // if (NaN_in_matrix (nrows, ncols, AbsResid, AbsResid_stride))
    //   cerr << "relative_residual: matrix AbsResid := A - Q*R contains a NaN" << endl;

    const magnitude_type absolute_residual = 
      local_frobenius_norm (nrows, ncols, &AbsResid[0], AbsResid_stride);

    // if (b_debug)
    //   {
    //     cerr << "In relative_residual:" << endl;
    //     cerr << "||Q||_2 = " << matrix_2norm(nrows, ncols, Q, ldq) << endl;
    //     cerr << "||R||_2 = " << matrix_2norm(ncols, ncols, R, ldr) << endl;
    //     cerr << "||A - QR||_2 = " << absolute_residual << endl;
    //   }

    return absolute_residual / A_norm_F;
  }

  /// Test accuracy of the computed QR factorization of the matrix A
  ///
  /// \param nrows [in] Number of rows in the A and Q matrices;
  ///   nrows >= ncols >= 1
  /// \param ncols [in] Number of columns in the A, Q, and R matrices;
  ///   nrows >= ncols >= 1
  /// \param A [in] Column-oriented nrows by ncols matrix with leading
  ///   dimension lda
  /// \param lda [in] Leading dimension of the matrix A; lda >= nrows
  /// \param Q [in] Column-oriented nrows by ncols matrix with leading
  ///   dimension ldq; computed Q factor of A
  /// \param ldq [in] Leading dimension of the matrix Q; ldq >= nrows
  /// \param R [in] Column-oriented upper triangular ncols by ncols 
  ///   matrix with leading dimension ldr; computed R factor of A
  /// \param ldr [in] Leading dimension of the matrix R; ldr >= ncols
  /// \return $\| A - Q R \|_F$, $\| I - Q^* Q \|_F$, and $\|A\|_F$.
  ///   The first is the residual of the QR factorization, the second 
  ///   a measure of the orthogonality of the resulting Q factor, and
  ///   the third an appropriate scaling factor if we want to compute 
  ///   the relative residual.  All are measured in the Frobenius 
  ///   (square root of (sum of squares of the matrix entries) norm.
  ///
  /// \note The reason for the elaborate "magnitude_type" construction
  /// is because this function returns norms, and norms always have
  /// real-valued type.  Scalar may be complex.  We could simply set
  /// the imaginary part to zero, but it seems more sensible to
  /// enforce the norm's value property in the type system.  Besides,
  /// one could imagine more elaborate Scalars (like rational
  /// functions, which do form a field) that have different plausible
  /// definitions of magnitude -- this is not just a problem for
  /// complex numbers (that are isomorphic to pairs of real numbers).
  template< class Ordinal, class Scalar >
  std::vector< typename ScalarTraits< Scalar >::magnitude_type >
  local_verify (const Ordinal nrows, 
		const Ordinal ncols, 
		const Scalar* const A, 
		const Ordinal lda,
		const Scalar* const Q, 
		const Ordinal ldq, 
		const Scalar* const R, 
		const Ordinal ldr)
  {
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
    std::vector< magnitude_type > results(3);
    // const bool A_contains_NaN = NaN_in_matrix (nrows, ncols, A, lda);
    // const bool Q_contains_NaN = NaN_in_matrix (nrows, ncols, Q, ldq);
    // const bool R_contains_NaN = NaN_in_matrix (ncols, ncols, R, ldr);

    results[0] = localResidual (nrows, ncols, A, lda, Q, ldq, R, ldr);
    results[1] = localOrthogonality (nrows, ncols, Q, ldq);
    results[2] = local_frobenius_norm (nrows, ncols, A, lda);

    return results;
  }

} // namespace TSQR

#endif // __TSQR_Tsqr_LocalVerify_hpp
