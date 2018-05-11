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

/// \file Tsqr_CombineNative.hpp
/// \brief Interface to C++ back end of \c TSQR::Combine.
///
#ifndef __TSQR_CombineNative_hpp
#define __TSQR_CombineNative_hpp

#include <Teuchos_LAPACK.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Tsqr_ApplyType.hpp"
#include "Tsqr_CombineDefault.hpp"

#include "Kokkos_Core.hpp"
#include "KokkosBlas2_gemv.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace TSQR {

  /// \class CombineNative
  /// \brief Interface to C++ back end of TSQR::Combine
  ///
  /// \c TSQR::Combine has three implementations: \c CombineDefault,
  /// CombineNative, and \c CombineFortran.  CombineNative,
  /// implemented in this file, is a fully C++ (therefore "native," as
  /// opposed to \c CombineFortran (implemented in Fortran) or \c
  /// CombineNative (implemented by wrappers around LAPACK calls))
  /// implementation.
  ///
  /// \warning CombineNative has no complex-arithmetic implementation
  ///   yet.  It's not hard to implement this (use LAPACK's ZGEQR2(P)
  ///   and ZUNM2R as models), but it will take time that the author
  ///   doesn't have at the moment.
  ///
  template< class Ordinal, class Scalar, bool isComplex = Teuchos::ScalarTraits< Scalar >::isComplex >
  class CombineNative
  {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef Ordinal ordinal_type;

  private:
    typedef Teuchos::LAPACK<ordinal_type, scalar_type> lapack_type;
    typedef CombineDefault<ordinal_type, scalar_type> combine_default_type;

  public:

    CombineNative () {}

    /// Whether or not the QR factorizations computed by methods of
    /// this class produce an R factor with all nonnegative diagonal
    /// entries.  It depends on LAPACK because this implementation
    /// invokes one of {LARFGP, LARFP, LARFG} in order to compute
    /// Householder reflectors; only LAPACK versions >= 3.2 have one
    /// of {LARFGP, LARFP}, which is necessary to ensure that the BETA
    /// output of the function is always nonnegative.
    static bool QR_produces_R_factor_with_nonnegative_diagonal() {
      return /* lapack_type::QR_produces_R_factor_with_nonnegative_diagonal() */ false &&
        combine_default_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const Ordinal nrows,
                  const Ordinal ncols,
                  Scalar A[],
                  const Ordinal lda,
                  Scalar tau[],
                  Scalar work[]) const
    {
      return default_.factor_first (nrows, ncols, A, lda, tau, work);
    }

    void
    apply_first (const ApplyType& applyType,
                 const Ordinal nrows,
                 const Ordinal ncols_C,
                 const Ordinal ncols_A,
                 const Scalar A[],
                 const Ordinal lda,
                 const Scalar tau[],
                 Scalar C[],
                 const Ordinal ldc,
                 Scalar work[]) const
    {
      return default_.apply_first (applyType, nrows, ncols_C, ncols_A,
                                   A, lda, tau, C, ldc, work);
    }

    void
    apply_inner (const ApplyType& applyType,
                 const Ordinal m,
                 const Ordinal ncols_C,
                 const Ordinal ncols_Q,
                 const Scalar A[],
                 const Ordinal lda,
                 const Scalar tau[],
                 Scalar C_top[],
                 const Ordinal ldc_top,
                 Scalar C_bot[],
                 const Ordinal ldc_bot,
                 Scalar work[]) const;

    void
    factor_inner (const Ordinal m,
                  const Ordinal n,
                  Scalar R[],
                  const Ordinal ldr,
                  Scalar A[],
                  const Ordinal lda,
                  Scalar tau[],
                  Scalar work[]) const;

    void
    factor_pair (const Ordinal n,
                 Scalar R_top[],
                 const Ordinal ldr_top,
                 Scalar R_bot[],
                 const Ordinal ldr_bot,
                 Scalar tau[],
                 Scalar work[]) const;

    void
    apply_pair (const ApplyType& applyType,
                const Ordinal ncols_C,
                const Ordinal ncols_Q,
                const Scalar R_bot[],
                const Ordinal ldr_bot,
                const Scalar tau[],
                Scalar C_top[],
                const Ordinal ldc_top,
                Scalar C_bot[],
                const Ordinal ldc_bot,
                Scalar work[]) const;

  private:
    mutable combine_default_type default_;
  };


  //! Specialization of CombineNative for the real-arithmetic case.
  template< class Ordinal, class Scalar >
  class CombineNative< Ordinal, Scalar, false >
  {
  private:
    using memory_space = Kokkos::HostSpace;
#ifdef KOKKOS_ENABLE_SERIAL
    using execution_space = Kokkos::Serial;
#else // NOT KOKKOS_ENABLE_SERIAL
    using execution_space = Kokkos::HostSpace::execution_space;
#endif // KOKKOS_ENABLE_SERIAL

  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef Ordinal ordinal_type;
    using device_type = Kokkos::Device<execution_space, memory_space>;

  private:
    typedef Teuchos::LAPACK<ordinal_type, scalar_type> lapack_type;
    typedef CombineDefault<ordinal_type, scalar_type> combine_default_type;

    void
    GER (const magnitude_type alpha,
         const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& x,
         const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& y,
         const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& A) const;

    void
    LARFG (const Ordinal n,
           scalar_type* const alpha,
           const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& x,
           scalar_type* const tau) const
    {
      constexpr Ordinal incx {1};
      lapack_type ().LARFG (n, alpha, x.data (), incx, tau);
    }

    magnitude_type
    LAPY2 (const scalar_type& x, const scalar_type& y) const
    {
      using KAT = Kokkos::ArithTraits<scalar_type>;
      if (KAT::isNan (x)) {
        return x;
      }
      else if (KAT::isNan (y)) {
        return y;
      }
      else {
        const magnitude_type xabs = KAT::abs (x);
        const magnitude_type yabs = KAT::abs (y);
        const scalar_type w = xabs >= yabs ? xabs : yabs; // max (xabs, yabs);
        const scalar_type z = xabs <= yabs ? xabs : yabs; // min (xabs, yabs);

        if (z == KAT::zero ()) {
          return w;
        }
        else {
          const scalar_type z_div_w = z / w;
          return w * KAT::sqrt (KAT::one () + z_div_w * z_div_w);
        }
      }
    }

    void
    GEMV (const char trans[],
          const scalar_type alpha,
          const Kokkos::View<const scalar_type**, Kokkos::LayoutLeft, device_type>& A,
          const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& x,
          const scalar_type beta,
          const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& y) const;

    void
    factor_pair (const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& R_top,
                 const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& R_bot,
                 const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& tau_view,
                 const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& work_view) const;

    void
    factor_inner (const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& R_view,
                  const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& A_view,
                  const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& tau_view,
                  const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& work_view) const;

    void
    apply_pair (const ApplyType& applyType,
                const Kokkos::View<const scalar_type**, Kokkos::LayoutLeft, device_type>& R_bot, // ncols_Q
                const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& tau_view,
                const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& C_top, // ncols_C
                const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& C_bot,
                const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& work_view) const;

    void
    apply_inner (const ApplyType& applyType,
                 const Kokkos::View<const scalar_type**, Kokkos::LayoutLeft, device_type>& A,
                 const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& tau,
                 const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& C_top,
                 const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& C_bot,
                 const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& work) const;

  public:
    CombineNative () = default;

    static bool QR_produces_R_factor_with_nonnegative_diagonal() {
      return /* lapack_type::QR_produces_R_factor_with_nonnegative_diagonal() */ false &&
        combine_default_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const Ordinal nrows,
                  const Ordinal ncols,
                  Scalar A[],
                  const Ordinal lda,
                  Scalar tau[],
                  Scalar work[]) const
    {
      return default_.factor_first (nrows, ncols, A, lda, tau, work);
    }

    void
    apply_first (const ApplyType& applyType,
                 const Ordinal nrows,
                 const Ordinal ncols_C,
                 const Ordinal ncols_A,
                 const Scalar A[],
                 const Ordinal lda,
                 const Scalar tau[],
                 Scalar C[],
                 const Ordinal ldc,
                 Scalar work[]) const
    {
      return default_.apply_first (applyType, nrows, ncols_C, ncols_A,
                                   A, lda, tau, C, ldc, work);
    }

    void
    apply_inner (const ApplyType& applyType,
                 const Ordinal m,
                 const Ordinal ncols_C,
                 const Ordinal ncols_Q,
                 const Scalar A[],
                 const Ordinal lda,
                 const Scalar tau[],
                 Scalar C_top[],
                 const Ordinal ldc_top,
                 Scalar C_bot[],
                 const Ordinal ldc_bot,
                 Scalar work[]) const;

    void
    factor_inner (const Ordinal m,
                  const Ordinal n,
                  Scalar R[],
                  const Ordinal ldr,
                  Scalar A[],
                  const Ordinal lda,
                  Scalar tau[],
                  Scalar work[]) const;

    void
    factor_pair (const Ordinal n,
                 Scalar R_top[],
                 const Ordinal ldr_top,
                 Scalar R_bot[],
                 const Ordinal ldr_bot,
                 Scalar tau[],
                 Scalar work[]) const;

    void
    apply_pair (const ApplyType& applyType,
                const Ordinal ncols_C,
                const Ordinal ncols_Q,
                const scalar_type R_bot[],
                const Ordinal ldr_bot,
                const scalar_type tau[],
                scalar_type C_top[],
                const Ordinal ldc_top,
                scalar_type C_bot[],
                const Ordinal ldc_bot,
                scalar_type work[]) const;

  private:
    mutable combine_default_type default_;
  };


  /// "Forward declaration" for the complex-arithmetic case.
  ///
  template< class Ordinal, class Scalar >
  class CombineNative< Ordinal, Scalar, true >
  {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef Ordinal ordinal_type;

  private:
    typedef Teuchos::LAPACK<ordinal_type, scalar_type> lapack_type;
    typedef CombineDefault<ordinal_type, scalar_type> combine_default_type;

  public:
    CombineNative () {}

    static bool QR_produces_R_factor_with_nonnegative_diagonal() {
      return /* lapack_type::QR_produces_R_factor_with_nonnegative_diagonal() */ false &&
        combine_default_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const Ordinal nrows,
                  const Ordinal ncols,
                  Scalar A[],
                  const Ordinal lda,
                  Scalar tau[],
                  Scalar work[]) const
    {
      return default_.factor_first (nrows, ncols, A, lda, tau, work);
    }

    void
    apply_first (const ApplyType& applyType,
                 const Ordinal nrows,
                 const Ordinal ncols_C,
                 const Ordinal ncols_A,
                 const Scalar A[],
                 const Ordinal lda,
                 const Scalar tau[],
                 Scalar C[],
                 const Ordinal ldc,
                 Scalar work[]) const
    {
      return default_.apply_first (applyType, nrows, ncols_C, ncols_A,
                                   A, lda, tau, C, ldc, work);
    }

    void
    apply_inner (const ApplyType& applyType,
                 const Ordinal m,
                 const Ordinal ncols_C,
                 const Ordinal ncols_Q,
                 const Scalar A[],
                 const Ordinal lda,
                 const Scalar tau[],
                 Scalar C_top[],
                 const Ordinal ldc_top,
                 Scalar C_bot[],
                 const Ordinal ldc_bot,
                 Scalar work[]) const
    {
      return default_.apply_inner (applyType, m, ncols_C, ncols_Q,
                                   A, lda, tau,
                                   C_top, ldc_top, C_bot, ldc_bot,
                                   work);
    }

    void
    factor_inner (const Ordinal m,
                  const Ordinal n,
                  Scalar R[],
                  const Ordinal ldr,
                  Scalar A[],
                  const Ordinal lda,
                  Scalar tau[],
                  Scalar work[]) const
    {
      return default_.factor_inner (m, n, R, ldr, A, lda, tau, work);
    }

    void
    factor_pair (const Ordinal n,
                 Scalar R_top[],
                 const Ordinal ldr_top,
                 Scalar R_bot[],
                 const Ordinal ldr_bot,
                 Scalar tau[],
                 Scalar work[]) const
    {
      return default_.factor_pair (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);
    }

    void
    apply_pair (const ApplyType& applyType,
                const Ordinal ncols_C,
                const Ordinal ncols_Q,
                const Scalar R_bot[],
                const Ordinal ldr_bot,
                const Scalar tau[],
                Scalar C_top[],
                const Ordinal ldc_top,
                Scalar C_bot[],
                const Ordinal ldc_bot,
                Scalar work[]) const
    {
      return default_.apply_pair (applyType, ncols_C, ncols_Q,
                                  R_bot, ldr_bot, tau,
                                  C_top, ldc_top, C_bot, ldc_bot,
                                  work);
    }

  private:
    mutable combine_default_type default_;
  };


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  GER (const magnitude_type alpha,
       const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& x,
       const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& y,
       const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& A) const
  {
    constexpr scalar_type ZERO {0.0};
    const Ordinal m = A.extent (0);
    const Ordinal n = A.extent (1);

    constexpr Ordinal incy {1};
    //Ordinal jy = (incy > 0) ? 1 : 1 - (n-1) * incy;
    Ordinal jy = 1;

    for (Ordinal j = 0; j < n; ++j) {
      if (y[jy-1] != ZERO) {
        const scalar_type temp = alpha * y[jy-1];
        for (Ordinal i = 0; i < m; ++i) {
          A(i,j) = A(i,j) + x[i] * temp;
        }
      }
      jy += incy;
    }
  }


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  GEMV (const char trans[],
        const scalar_type alpha,
        const Kokkos::View<const scalar_type**, Kokkos::LayoutLeft, device_type>& A,
        const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& x,
        const scalar_type beta,
        const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& y) const
  {
    using y_vec_type = Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>;
    using x_vec_type = Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>;
    using range_type = std::pair<Ordinal, Ordinal>;

    const Ordinal m = A.extent (0);
    const Ordinal n = A.extent (1);

    const bool no_trans = (trans[0] == 'N' || trans[0] == 'n');
    x_vec_type x_view = Kokkos::subview (x, range_type (0, no_trans ? n : m));
    y_vec_type y_view = Kokkos::subview (y, range_type (0, no_trans ? m : n));

    KokkosBlas::gemv (trans, alpha, A, x_view, beta, y_view);
  }

  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  factor_inner (const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& R_view,
                const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& A_view,
                const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& tau_view,
                const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& work_view) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using range_type = std::pair<Ordinal, Ordinal>;
    constexpr scalar_type ZERO {0.0};
    constexpr scalar_type ONE {1.0};

    const Ordinal m = A_view.extent (0);
    const Ordinal n = A_view.extent (1);

    for (Ordinal k = 0; k < n; ++k) {
      work_view(k) = ZERO;
    }

    for (Ordinal k = 0; k < n-1; ++k) {
      Scalar& R_kk = R_view(k, k);
      auto A_1k = subview (A_view, ALL (), k);
      auto A_1kp1 = subview (A_view, range_type (0, m), range_type (k+1, n));

      this->LARFG (m + 1, &R_kk, A_1k, &tau_view[k]);
      this->GEMV ("T", ONE, A_1kp1, A_1k, ZERO, work_view);

      for (Ordinal j = k+1; j < n; ++j) {
        Scalar& R_kj = R_view(k, j);

        work_view(j-k-1) += R_kj;
        R_kj -= tau_view[k] * work_view(j-k-1);
      }
      this->GER (-tau_view[k], A_1k, work_view, A_1kp1);
    }
    Scalar& R_nn = R_view(n-1, n-1);
    auto A_1n = subview (A_view, ALL (), n-1);

    this->LARFG (m+1, &R_nn, A_1n, &tau_view[n-1]);
  }


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  factor_inner (const Ordinal m,
                const Ordinal n,
                Scalar R[],
                const Ordinal ldr,
                Scalar A[],
                const Ordinal lda,
                Scalar tau[],
                Scalar work[]) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using mat_type =
      Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>;
    using nonconst_vec_type =
      Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>;
    using range_type = std::pair<Ordinal, Ordinal>;

    mat_type A_full (A, lda, n);
    mat_type A_view = subview (A_full, range_type (0, m), ALL ());
    mat_type R_full (R, ldr, n);
    mat_type R_view = subview (R_full, range_type (0, n), ALL ());
    nonconst_vec_type tau_view (tau, n);
    nonconst_vec_type work_view (work, n);

    this->factor_inner (R_view, A_view, tau_view, work_view);
  }

  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  apply_inner (const ApplyType& applyType,
               const Kokkos::View<const scalar_type**, Kokkos::LayoutLeft, device_type>& A,
               const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& tau,
               const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& C_top,
               const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& C_bot,
               const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& work) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using const_vec_type =
      Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>;
    constexpr scalar_type ZERO {0.0};

    const Ordinal m = A.extent (0);
    const Ordinal ncols_Q = A.extent (1);
    const Ordinal ncols_C = C_top.extent (1);

    for (Ordinal i = 0; i < ncols_C; ++i) {
      work(i) = ZERO;
    }

    Ordinal j_start, j_end, j_step;
    if (applyType == ApplyType::NoTranspose) {
      j_start = ncols_Q - 1;
      j_end = -1; // exclusive
      j_step = -1;
    }
    else {
      j_start = 0;
      j_end = ncols_Q; // exclusive
      j_step = +1;
    }
    for (Ordinal j = j_start; j != j_end; j += j_step) {
      const_vec_type A_1j = subview (A, ALL (), j);

      //blas.GEMV ("T", m, ncols_C, ONE, C_bot, ldc_bot, A_1j, 1, ZERO, &y[0], 1);
      for (Ordinal i = 0; i < ncols_C; ++i) {
        work(i) = ZERO;
        for (Ordinal k = 0; k < m; ++k) {
          work(i) += A_1j(k) * C_bot(k, i);
        }
        work(i) += C_top(j, i);
      }
      for (Ordinal k = 0; k < ncols_C; ++k) {
        C_top(j, k) -= tau[j] * work(k);
      }

      this->GER (-tau[j], A_1j, work, C_bot);
    }
  }

  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  apply_inner (const ApplyType& applyType,
               const Ordinal m,
               const Ordinal ncols_C,
               const Ordinal ncols_Q,
               const Scalar A[],
               const Ordinal lda,
               const Scalar tau[],
               Scalar C_top[],
               const Ordinal ldc_top,
               Scalar C_bot[],
               const Ordinal ldc_bot,
               Scalar work[]) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using const_mat_type =
      Kokkos::View<const scalar_type**, Kokkos::LayoutLeft, device_type>;
    using nonconst_mat_type =
      Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>;
    using const_vec_type =
      Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>;
    using nonconst_vec_type =
      Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>;
    using range_type = std::pair<Ordinal, Ordinal>;

    const_mat_type A_full (A, lda, ncols_Q);
    auto A_view = subview (A_full, range_type (0, m), ALL ());
    nonconst_mat_type C_top_full (C_top, ldc_top, ncols_C);
    auto C_top_view = subview (C_top_full, range_type (0, m), ALL ());
    nonconst_mat_type C_bot_full (C_bot, ldc_bot, ncols_C);
    auto C_bot_view = subview (C_bot_full, range_type (0, m), ALL ());
    const_vec_type tau_view (tau, ncols_Q);
    nonconst_vec_type work_view (work, ncols_C);

    this->apply_inner (applyType, A_view, tau_view, C_top_view, C_bot_view, work_view);
  }


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  factor_pair (const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& R_top,
               const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& R_bot,
               const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& tau_view,
               const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& work_view) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using range_type = std::pair<Ordinal, Ordinal>;
    constexpr scalar_type ZERO {0.0};
    constexpr scalar_type ONE {1.0};

    const Ordinal n = R_top.extent (0);
    for (Ordinal k = 0; k < n; ++k) {
      work_view(k) = ZERO;
    }

    for (Ordinal k = 0; k < n-1; ++k) {
      scalar_type& R_top_kk = R_top(k, k);
      auto R_bot_1k = subview (R_bot, ALL (), k);
      auto R_bot_1kp1 = subview (R_bot, range_type (0, k+1), range_type (k+1, n));

      // k+2: 1 element in R_top (R_top(k,k)), and k+1 elements in
      // R_bot (R_bot(1:k,k), in 1-based indexing notation).
      this->LARFG (k+2, &R_top_kk, R_bot_1k, &tau_view[k]);
      // One-based indexing, Matlab version of the GEMV call below:
      // work(1:k) := R_bot(1:k,k+1:n)' * R_bot(1:k,k)

      this->GEMV ("T", ONE, R_bot_1kp1, R_bot_1k, ZERO, work_view);

      for (Ordinal j = k+1; j < n; ++j) {
        scalar_type& R_top_kj = R_top(k, j);
        work_view(j-k-1) += R_top_kj;
        R_top_kj -= tau_view[k] * work_view(j-k-1);
      }
      this->GER (-tau_view[k], R_bot_1k, work_view, R_bot_1kp1);
    }
    scalar_type& R_top_nn = R_top(n-1, n-1);
    auto R_bot_1n = subview (R_bot, ALL (), n-1);

    // n+1: 1 element in R_top (n,n), and n elements in R_bot (the
    // whole last column).
    this->LARFG (n+1, &R_top_nn, R_bot_1n, &tau_view[n-1]);
  }


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  factor_pair (const Ordinal n,
               Scalar R_top[],
               const Ordinal ldr_top,
               Scalar R_bot[],
               const Ordinal ldr_bot,
               Scalar tau[],
               Scalar work[]) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using range_type = std::pair<Ordinal, Ordinal>;

    Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type> R_top_full (R_top, ldr_top, n);
    Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type> R_bot_full (R_bot, ldr_bot, n);
    Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type> tau_view (tau, n);
    Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type> work_view (work, n);

    if (ldr_top == n) {
      if (ldr_bot == n) {
        this->factor_pair (R_top_full, R_bot_full, tau_view, work_view);
      }
      else {
        auto R_bot_view = subview (R_bot_full, range_type (0, n), ALL ());
        this->factor_pair (R_top_full, R_bot_view, tau_view, work_view);
      }
    }
    else {
      auto R_top_view = subview (R_top_full, range_type (0, n), ALL ());
      if (ldr_bot == n) {
        this->factor_pair (R_top_view, R_bot_full, tau_view, work_view);
      }
      else {
        auto R_bot_view = subview (R_bot_full, range_type (0, n), ALL ());
        this->factor_pair (R_top_view, R_bot_view, tau_view, work_view);
      }
    }
  }


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  apply_pair (const ApplyType& applyType,
              const Ordinal ncols_C,
              const Ordinal ncols_Q,
              const scalar_type R_bot[],
              const Ordinal ldr_bot,
              const scalar_type tau[],
              scalar_type C_top[],
              const Ordinal ldc_top,
              scalar_type C_bot[],
              const Ordinal ldc_bot,
              scalar_type work[]) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using range_type = std::pair<Ordinal, Ordinal>;
    using const_mat_type =
      Kokkos::View<const scalar_type**, Kokkos::LayoutLeft, device_type>;
    using nonconst_mat_type =
      Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>;
    using const_vec_type =
      Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>;
    using nonconst_vec_type =
      Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>;

    const_mat_type R_bot_full (R_bot, ldr_bot, ncols_Q);
    nonconst_mat_type C_top_full (C_top, ldc_top, ncols_C);
    nonconst_mat_type C_bot_full (C_bot, ldc_bot, ncols_C);
    const_vec_type tau_view (tau, ncols_Q);
    nonconst_vec_type work_view (work, ncols_C);

    auto R_bot_view = subview (R_bot_full, range_type (0, ncols_Q), ALL ());
    auto C_top_view = subview (C_top_full, range_type (0, ncols_C), ALL ());
    auto C_bot_view = subview (C_bot_full, range_type (0, ncols_C), ALL ());
    this->apply_pair (applyType, R_bot_view, tau_view, C_top_view, C_bot_view, work_view);
  }

  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  apply_pair (const ApplyType& applyType,
              const Kokkos::View<const scalar_type**, Kokkos::LayoutLeft, device_type>& R_bot, // ncols_Q
              const Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>& tau_view,
              const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& C_top, // ncols_C
              const Kokkos::View<scalar_type**, Kokkos::LayoutLeft, device_type>& C_bot,
              const Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type>& work_view) const
  {
    using const_vec_type =
      Kokkos::View<const scalar_type*, Kokkos::LayoutLeft, device_type>;
    constexpr scalar_type ZERO {0.0};
    const Ordinal ncols_C = C_top.extent (1);
    const Ordinal ncols_Q = R_bot.extent (1);

    Ordinal j_start, j_end, j_step;
    if (applyType == ApplyType::NoTranspose) {
      j_start = ncols_Q - 1;
      j_end = -1; // exclusive
      j_step = -1;
    }
    else {
      j_start = 0;
      j_end = ncols_Q; // exclusive
      j_step = +1;
    }
    for (Ordinal j_Q = j_start; j_Q != j_end; j_Q += j_step) {
      // Using Householder reflector stored in column j_Q of R_bot
      const_vec_type R_bot_col = Kokkos::subview (R_bot, Kokkos::ALL (), j_Q);

      // In 1-based indexing notation, with k in 1, 2, ..., ncols_C
      // (inclusive): (Output is length ncols_C row vector)
      //
      // work(1:j) := R_bot(1:j,j)' * C_bot(1:j, 1:ncols_C) - C_top(j, 1:ncols_C)
      for (Ordinal j_C = 0; j_C < ncols_C; ++j_C) {
        // For each column j_C of [C_top; C_bot], update row j_Q
        // of C_top and rows 1:j_Q of C_bot.  (Again, this is in
        // 1-based indexing notation.

        scalar_type work_j_C = ZERO;
        const_vec_type C_bot_col = Kokkos::subview (C_bot, Kokkos::ALL (), j_C);

        for (Ordinal k = 0; k <= j_Q; ++k)
          work_j_C += R_bot_col(k) * C_bot_col(k);

        work_j_C += C_top(j_Q, j_C);
        work_view(j_C) = work_j_C;
      }
      for (Ordinal j_C = 0; j_C < ncols_C; ++j_C) {
        C_top(j_Q, j_C) -= tau_view[j_Q] * work_view(j_C);
      }
      this->GER (-tau_view[j_Q], R_bot_col, work_view, C_bot);
    }
  }
} // namespace TSQR



#endif // __TSQR_CombineNative_hpp
