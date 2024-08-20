// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Tsqr_LocalVerify_hpp
#define __TSQR_Tsqr_LocalVerify_hpp

#include "Tsqr_Util.hpp"
#include "Tsqr_Impl_SystemBlas.hpp"
#include "Tsqr_Matrix.hpp"
#include <cmath>
#include <limits>
#include <utility> // std::pair, std::make_pair
#include <vector>

namespace TSQR {
  template<class Ordinal, class Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  local_frobenius_norm (const Ordinal nrows_local,
                        const Ordinal ncols,
                        const Scalar A_local[],
                        const Ordinal lda_local)
  {
    using STS = Teuchos::ScalarTraits<Scalar>;
    using magnitude_type = typename STS::magnitudeType;

    // FIXME (mfh 22 Apr 2010) This function does no scaling of
    // intermediate quantities, so it might overflow unnecessarily.
    magnitude_type result (0);
    for (Ordinal j = 0; j < ncols; ++j) {
      const Scalar* const cur_col = &A_local[j*lda_local];
      for (Ordinal i = 0; i < nrows_local; ++i) {
        const auto abs_xi = STS::magnitude (cur_col[i]);
        result = result + abs_xi * abs_xi;
      }
    }
    // FIXME (mfh 14 Oct 2014) Should we use std::sqrt or even
    // STS::squareroot here instead?
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
#ifdef __CUDACC__
              if (isnan (A[i + j*lda]))
#else
              if (std::isnan (A[i + j*lda]))
#endif
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
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  localOrthogonality (const Ordinal nrows,
                      const Ordinal ncols,
                      const Scalar Q[],
                      const Ordinal ldq)
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    const Scalar ZERO {};
    const Scalar ONE (1.0);

    Impl::SystemBlas<Scalar> blas;

    std::vector<Scalar> AbsOrthog (ncols * ncols, std::numeric_limits<Scalar>::quiet_NaN());
    const Ordinal AbsOrthog_stride = ncols;

    // Compute AbsOrthog := Q' * Q - I.  First, compute Q' * Q:
    if (STS::isComplex) {
      blas.GEMM (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, ncols, ncols, nrows,
                 ONE, Q, ldq, Q, ldq, ZERO, &AbsOrthog[0], AbsOrthog_stride);
    }
    else {
      blas.GEMM (Teuchos::TRANS, Teuchos::NO_TRANS, ncols, ncols, nrows,
                 ONE, Q, ldq, Q, ldq, ZERO, &AbsOrthog[0], AbsOrthog_stride);
    }

    // Now, compute (Q^T*Q) - I.
    for (Ordinal j = 0; j < ncols; ++j) {
      AbsOrthog[j + j*AbsOrthog_stride] = AbsOrthog[j + j*AbsOrthog_stride] - ONE;
    }

    // Now AbsOrthog == Q^T * Q - I.  Compute and return its Frobenius norm.
    return local_frobenius_norm (ncols, ncols, &AbsOrthog[0], AbsOrthog_stride);
  }



  template< class Ordinal, class Scalar >
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  local_relative_orthogonality (const Ordinal nrows,
                                const Ordinal ncols,
                                const Scalar Q[],
                                const Ordinal ldq,
                                const typename Teuchos::ScalarTraits<Scalar>::magnitudeType A_norm_F)
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;
    const Scalar ZERO {};
    const Scalar ONE (1.0);

    const bool relative = false; // whether to scale $\|I-Q^T*Q\|_F$ by $\|A\|_F$
    Impl::SystemBlas<Scalar> blas;

    std::vector<Scalar> AbsOrthog (ncols * ncols, std::numeric_limits<Scalar>::quiet_NaN());
    const Ordinal AbsOrthog_stride = ncols;

    // Compute AbsOrthog := Q' * Q - I.  First, compute Q' * Q:
    if (STS::isComplex) {
      blas.GEMM (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, ncols, ncols, nrows,
                 ONE, Q, ldq, Q, ldq, ZERO, &AbsOrthog[0], AbsOrthog_stride);
    }
    else {
      blas.GEMM (Teuchos::TRANS, Teuchos::NO_TRANS, ncols, ncols, nrows,
                 ONE, Q, ldq, Q, ldq, ZERO, &AbsOrthog[0], AbsOrthog_stride);
    }

    // Now, compute (Q^T*Q) - I.
    for (Ordinal j = 0; j < ncols; ++j) {
      AbsOrthog[j + j*AbsOrthog_stride] = AbsOrthog[j + j*AbsOrthog_stride] - ONE;
    }

    // Now AbsOrthog == Q^T * Q - I.  Compute its Frobenius norm.
    const magnitude_type AbsOrthog_norm_F =
      local_frobenius_norm (ncols, ncols, &AbsOrthog[0], AbsOrthog_stride);

    // Return the orthogonality measure
    return relative ? (AbsOrthog_norm_F / A_norm_F) : AbsOrthog_norm_F;
  }


  template< class Ordinal, class Scalar >
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  localResidual (const Ordinal nrows,
                 const Ordinal ncols,
                 const Scalar A[],
                 const Ordinal lda,
                 const Scalar Q[],
                 const Ordinal ldq,
                 const Scalar R[],
                 const Ordinal ldr)
  {
    using Teuchos::NO_TRANS;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    MatView<Ordinal, const Scalar> A_view (nrows, ncols, A, lda);
    Matrix<Ordinal, Scalar> AbsResid (nrows, ncols,
      std::numeric_limits<Scalar>::quiet_NaN ());
    Impl::SystemBlas<Scalar> blas;
    const magnitude_type ONE (1);

    // A_copy := A_copy - Q * R
    deep_copy (AbsResid, A_view);
    blas.GEMM (NO_TRANS, NO_TRANS, nrows, ncols, ncols, -ONE, Q, ldq, R, ldr,
               ONE, AbsResid.data(), AbsResid.stride(1));

    return local_frobenius_norm (nrows, ncols, AbsResid.data(),
                                 AbsResid.stride(1));
  }


  template< class Ordinal, class Scalar >
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  local_relative_residual (const Ordinal nrows,
                           const Ordinal ncols,
                           const Scalar A[],
                           const Ordinal lda,
                           const Scalar Q[],
                           const Ordinal ldq,
                           const Scalar R[],
                           const Ordinal ldr,
                           const typename Teuchos::ScalarTraits<Scalar>::magnitudeType A_norm_F)
  {
    using Teuchos::NO_TRANS;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    MatView<Ordinal, const Scalar> A_view (nrows, ncols, A, lda);
    Matrix<Ordinal, Scalar> AbsResid
      (nrows, ncols, std::numeric_limits<Scalar>::quiet_NaN ());
    deep_copy (AbsResid, A);

    // A_copy := A_copy - Q * R
    Impl::SystemBlas<Scalar> blas;
    const magnitude_type ONE (1.0);
    blas.GEMM (NO_TRANS, NO_TRANS, nrows, ncols, ncols,
               -ONE, Q, ldq, R, ldr,
               ONE, AbsResid.data(), AbsResid.stride(1));

    const magnitude_type absolute_residual =
      local_frobenius_norm (nrows, ncols, AbsResid.data(),
                            AbsResid.stride(1));
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
  std::vector< typename Teuchos::ScalarTraits<Scalar>::magnitudeType >
  local_verify (const Ordinal nrows,
                const Ordinal ncols,
                const Scalar* const A,
                const Ordinal lda,
                const Scalar* const Q,
                const Ordinal ldq,
                const Scalar* const R,
                const Ordinal ldr)
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;
    std::vector<magnitude_type> results (3);
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
