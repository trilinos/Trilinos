// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Tsqr_GlobalVerify_hpp
#define __TSQR_Tsqr_GlobalVerify_hpp

#include "Tsqr_LocalVerify.hpp"
#include "Tsqr_MessengerBase.hpp"
#include "Tsqr_Util.hpp"
#include "Tsqr_Impl_SystemBlas.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <utility> // std::pair
#include <vector>

namespace TSQR {
  /// \class GlobalSummer
  ///
  /// Compute a global sum of (magnitudes of) Scalar values, returning
  /// a magnitude_type.
  ///
  /// \note Unfortunately, you need C++11 support to have default
  ///   template arguments of template functions.  Otherwise we would
  ///   make this a template function and set the default value of
  ///   isComplex to Teuchos::ScalarTraits<Scalar>::isComplex.  Also,
  ///   C++ (before C++11) doesn't like partial specialization of
  ///   template functions.  So, we had to make this a class.
  template<class Scalar, bool isComplex = Teuchos::ScalarTraits< Scalar >::isComplex>
  class GlobalSummer {
  public:
    typedef Scalar scalar_type;
    typedef Teuchos::ScalarTraits< Scalar > STS;
    typedef typename STS::magnitudeType magnitude_type;

    static magnitude_type
    sum (const Scalar& localSum,
         MessengerBase<Scalar>* const messenger);
  };

  // Complex-arithmetic "forward declaration"
  template<class Scalar>
  class GlobalSummer<Scalar, true> {
  public:
    typedef Scalar scalar_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    static magnitude_type
    sum (const Scalar& localSum,
         MessengerBase<Scalar>* const messenger);
  };

  // Real-arithmetic "forward declaration"
  template<class Scalar>
  class GlobalSummer<Scalar, false> {
  public:
    typedef Scalar scalar_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    static magnitude_type
    sum (const Scalar& localSum,
         MessengerBase<Scalar>* const messenger);
  };

  // Complex-arithmetic case
  template<class Scalar>
  typename GlobalSummer<Scalar, true>::magnitude_type
  GlobalSummer<Scalar, true>::sum (const Scalar& localSum,
                                   MessengerBase<Scalar>* const messenger)
  {
    // In order to use a MessengerBase<Scalar> on magnitude_type
    // values, we have to convert local_result to a Scalar, and then
    // convert back the result.  We convert by setting the real
    // component of the Scalar to the magnitude_type.  This isn't
    // guaranteed to work if magnitude_type has a greater dynamic
    // range than Scalar.  That's possible, but that's not how we do
    // things with ScalarTraits< std::complex< T > >, and that's not
    // how LAPACK does it either, so it's fair to assume that
    // magnitude_type and the individual components of Scalar have the
    // same dynamic range.
    const magnitude_type localSumAbs = STS::magnitude (localSum);
    const Scalar localSumAsScalar (localSumAbs, magnitude_type(0));
    const Scalar globalSumAsScalar = messenger->globalSum (localSumAsScalar);
    const magnitude_type globalSum = STS::magnitude (globalSumAsScalar);
    return globalSum;
  }

  // Real-arithmetic case
  template<class Scalar>
  typename GlobalSummer<Scalar, false>::magnitude_type
  GlobalSummer<Scalar, false>::sum (const Scalar& localSum,
                                    MessengerBase<Scalar>* const messenger)
  {
    const Scalar localSumAsScalar (localSum);
    const Scalar globalSumAsScalar = messenger->globalSum (localSumAsScalar);
    const magnitude_type globalSum = STS::magnitude (globalSumAsScalar);
    return globalSum;
  }

  template<class LocalOrdinal, class Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  global_frobenius_norm (const LocalOrdinal nrows_local,
                         const LocalOrdinal ncols,
                         const Scalar A_local[],
                         const LocalOrdinal lda_local,
                         MessengerBase<Scalar>* const messenger)
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    magnitude_type localResult {};
    for (LocalOrdinal j = 0; j < ncols; j++) {
      const Scalar* const cur_col = &A_local[j*lda_local];
      for (LocalOrdinal i = 0; i < nrows_local; ++i) {
        const magnitude_type abs_xi = STS::magnitude (cur_col[i]);
        localResult = localResult + abs_xi * abs_xi;
      }
    }
    // GlobalSummmer() is a hack to let us use a Scalar - type
    // MessengerBase with magnitude_type inputs and outputs.
    // Otherwise we would need to carry around a
    // MessengerBase<magnitude_type> object as well.
    const magnitude_type globalResult =
      GlobalSummer<Scalar, STS::isComplex>::sum (localResult, messenger);
    return sqrt (globalResult);
  }

  template<class LocalOrdinal, class Scalar>
  std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
  global_verify (const LocalOrdinal nrows_local,
                 const LocalOrdinal ncols,
                 const Scalar A_local[],
                 const LocalOrdinal lda_local,
                 const Scalar Q_local[],
                 const LocalOrdinal ldq_local,
                 const Scalar R[],
                 const LocalOrdinal ldr,
                 MessengerBase<Scalar>* const messenger)
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;
    using Teuchos::CONJ_TRANS;
    using Teuchos::NO_TRANS;
    using Teuchos::TRANS;
    using std::make_pair;
    using std::pair;
    using std::vector;

    const magnitude_type ZERO {};
    const magnitude_type ONE (1.0);
    Impl::SystemBlas<Scalar> blas;

    //
    // Compute $\| I - Q^T * Q \|_F$
    //

    // Compute Q_local^T * Q_local (this node's component of Q^T*Q)
    vector<Scalar> Temp (ncols*ncols, STS::nan());
    const LocalOrdinal ld_temp = ncols;

    if (STS::isComplex)
      blas.GEMM (CONJ_TRANS, NO_TRANS, ncols, ncols, nrows_local,
                 ONE, Q_local, ldq_local, Q_local, ldq_local,
                 ZERO, &Temp[0], ld_temp);
    else
      blas.GEMM (TRANS, NO_TRANS, ncols, ncols, nrows_local,
                 ONE, Q_local, ldq_local, Q_local, ldq_local,
                 ZERO, &Temp[0], ld_temp);

    // Reduce over all the processors to get the global Q^T*Q in Temp2.
    vector<Scalar> Temp2 (ncols*ncols, STS::nan());
    messenger->globalVectorSum (&Temp[0], &Temp2[0], ncols*ncols);

    // Compute I-(Q^T*Q) redundantly on all processors
    for (LocalOrdinal j = 0; j < ncols; j++)
      Temp2[j + j*ld_temp] = ONE - Temp2[j + j*ld_temp];

    // Compute the Frobenius norm of I - Q^T*Q, redundantly on all processors.
    const magnitude_type Orthog_F =
      local_frobenius_norm (ncols, ncols, &Temp2[0], ld_temp);

    // Compute the Frobenius norm of A.
    const magnitude_type A_F =
      global_frobenius_norm (nrows_local, ncols, &A_local[0], lda_local, messenger);

    //
    // Compute $\| A - Q*R \|_F$
    //

    vector<Scalar> Resid (nrows_local * ncols, STS::nan());
    const LocalOrdinal ld_resid = nrows_local;
    MatView<LocalOrdinal, Scalar> Resid_view
      (nrows_local, ncols, Resid.data (), ld_resid);
    MatView<LocalOrdinal, const Scalar> A_view
      (nrows_local, ncols, A_local, lda_local);
    deep_copy (Resid_view, A_view);

    // Resid := Resid - Q*R
    blas.GEMM (NO_TRANS, NO_TRANS, nrows_local, ncols, ncols,
               -ONE, Q_local, ldq_local, R, ldr,
               ONE, Resid.data(), ld_resid);

    const magnitude_type Resid_F =
      global_frobenius_norm (nrows_local, ncols, Resid.data(),
                             ld_resid, messenger);
    vector<magnitude_type> results (3);
    results[0] = Resid_F;
    results[1] = Orthog_F;
    results[2] = A_F;
    return results;
  }

} // namespace TSQR

#endif // __TSQR_Tsqr_GlobalVerify_hpp

