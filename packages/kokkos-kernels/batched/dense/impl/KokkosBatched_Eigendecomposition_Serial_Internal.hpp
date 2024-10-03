//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef __KOKKOSBATCHED_EIGENDECOMPOSITION_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_EIGENDECOMPOSITION_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_SetIdentity_Internal.hpp"
#include "KokkosBatched_SetTriangular_Internal.hpp"
#include "KokkosBatched_Normalize_Internal.hpp"
#include "KokkosBatched_Hessenberg_Serial_Internal.hpp"
#include "KokkosBatched_ApplyQ_Serial_Internal.hpp"
#include "KokkosBatched_Schur_Serial_Internal.hpp"
#include "KokkosBatched_RightEigenvectorFromSchur_Serial_Internal.hpp"
#include "KokkosBatched_LeftEigenvectorFromSchur_Serial_Internal.hpp"
#include "KokkosBatched_Gemm_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================

struct SerialEigendecompositionInternal {
  /// Given a general nonsymmetric matrix A (m x m), it performs
  /// eigendecomposition of the matrix.
  ///
  /// Parameters:
  ///   [in]m
  ///     A dimension of the square matrix H.
  ///   [in/out]A, [in]as0, [in]as1
  ///     Real general nonsymmetric matrix A(m x m) with strides as0 and as1.
  ///     A is first condensed to a upper Hessenberg form. Then, the Francis
  ///     double shift QR algorithm is applied to compute its Schur form.
  ///     On exit, A stores a quasi upper triangular matrix of the Schur
  ///     decomposition.
  ///   [out]er, [in]ers, [out]ei, [in]eis
  ///     A complex vector er(m)+ei(m)i with a stride ers and eis to store
  ///     computed eigenvalues. For a complex eigen pair, it stores a+bi and
  ///     a-bi consecutively.
  ///   [out]UL, [in]uls0, [in] uls1
  ///     Left eigenvectors UL(m x m) with strides uls0 and uls1. When UL is
  ///     NULL, the left eigenvectors are not computed.
  ///   [out]UR, [in]urs0, [in] urs1
  ///     Right eigenvectors UR(m x m) with strides urs0 and urs1. When UR is
  ///     NULL, the right eigenvectors are not computed.
  ///   [out]w, [in]wlen
  ///     Workspace
  template <typename RealType>
  KOKKOS_INLINE_FUNCTION static int device_invoke(const int m, RealType* A, const int as0, const int as1, RealType* er,
                                                  const int ers, RealType* ei, const int eis, RealType* UL,
                                                  const int uls0, const int uls1, RealType* UR, const int urs0,
                                                  const int urs1, RealType* w, const int wlen) {
    /// until debugging is done, comment out the code
    /// testing happens only for TPLs on host.
    static_assert(false,
                  "Serial eigendecomposition on device and/or without LAPACK "
                  "is not implemented yet");
    //       typedef RealType real_type;
    //       typedef Kokkos::ArithTraits<real_type> ats;

    //       const real_type one(1), zero(0), tol = 1e2*ats::epsilon();
    //       //const Kokkos::pair<real_type,real_type> identity(one, zero);

    //       /// step 0: input checking
    //       assert( (wlen >= (2*m*m+5*m)) && "Eigendecomposition: workspace
    //       size is too small"); real_type *w_now = w; int wlen_now = wlen;
    //       assert( (wlen_now >= 0) && "Eigendecomposition: workspace size is
    //       negative");

    //       const bool is_UL = UL != NULL, is_UR = UR != NULL;
    //       assert((is_UL || is_UR) && "Eigendecomposition: neither left nor
    //       right eigenvectors were requested. Use SerialEigenvalueInternal
    //       instead.");

    //       real_type *QZ = w_now; w_now += (m*m); wlen_now -= (m*m);
    //       assert( (wlen_now >= 0) && "Eigendecomposition: QZ allocation
    //       fails");

    //       const int qzs0 = m, qzs1 = 1, qzs = m+1; /// row major
    //       const int as = as0+as1;

    //       /// step 1: Hessenberg reduction A = Q H Q^H
    //       ///         Q is stored in QZ
    //
    ////////////////////////////////////////////////////////////////////////////
    // DO NOT USE
    //
    //     #ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    //     <host code>
    //     #else
    //     <device code>
    //     #endif
    //
    // DO THIS INSTEAD
    //
    //     KOKKOS_IF_ON_HOST((<host code>))
    //     KOKKOS_IF_ON_DEVICE((<device code>))
    //
    ////////////////////////////////////////////////////////////////////////////
    // #if (defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && (__INTEL_MKL__ >= 2018)) &&
    // defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    //       {
    //         real_type *t  = w_now; w_now += m; wlen_now -= m;
    //         assert( (wlen_now >= 0) && "Eigendecomposition: Hessenberg
    //         reduction workspace t allocation fail");

    //         if (as0 == 1 || as1 == 1) { /// if mkl can be interfaced, use it
    //           const auto matrix_layout = ( as1 == 1 ? LAPACK_ROW_MAJOR :
    //           LAPACK_COL_MAJOR ); LAPACKE_dgehrd(matrix_layout, m, 1, m, A,
    //           m, t);

    //           SerialSetIdentityInternal::invoke(m, QZ, qzs0, qzs1);
    //           LAPACKE_dorghr(matrix_layout, m, 1, m, QZ, m, t);
    //         } else { /// for arbitrary strides, there is no choice to use
    //         tpls
    //           real_type *ww = w_now; w_now += m; wlen_now -= m;
    //           assert( (wlen_now >= 0) && "Eigendecomposition: Hessenberg
    //           reduction workspace ww allocation fail");

    //           SerialHessenbergInternal::invoke(m, m,
    //                                            A, as0, as1,
    //                                            t, 1,
    //                                            ww);

    //           SerialSetIdentityInternal::invoke(m, QZ, qzs0, qzs1);
    //           SerialApplyQ_LeftForwardInternal::invoke(m-1, m-1, m-1,
    //                                                    A+as0, as0, as1,
    //                                                    t, 1,
    //                                                    QZ+qzs, qzs0, qzs1,
    //                                                    ww);
    //           /// recovery of workspace for ww
    //           w_now -= m; wlen_now += m;
    //         }
    //         /// recovery of workspace for t
    //         w_now -= m; wlen_now += m;

    //         /// clean up H
    //         SerialSetLowerTriangularInternal::invoke(m, m,
    //                                                  2,
    //                                                  zero,
    //                                                  A, as0, as1);
    //       }
    // #else
    //       {
    //         real_type *t  = w_now; w_now += m; wlen_now -= m;
    //         real_type *ww = w_now; w_now += m; wlen_now -= m;
    //         assert( (wlen_now >= 0) && "Eigendecomposition: Hessenberg
    //         reduction workspace t and ww allocation fail");

    //         SerialHessenbergInternal::invoke(m, m,
    //                                          A, as0, as1,
    //                                          t, 1,
    //                                          ww);

    //         SerialSetIdentityInternal::invoke(m, QZ, qzs0, qzs1);
    //         SerialApplyQ_LeftForwardInternal::invoke(m-1, m-1, m-1,
    //                                                  A+as0, as0, as1,
    //                                                  t, 1,
    //                                                  QZ+qzs, qzs0, qzs1,
    //                                                  ww);

    //         /// clean up H
    //         SerialSetLowerTriangularInternal::invoke(m, m,
    //                                                  2,
    //                                                  zero,
    //                                                  A, as0, as1);

    //         /// recover workspace
    //         w_now -= (2*m); wlen_now += (2*m);
    //       }
    // #endif
    //       /// step 2: Schur decomposition H = Z T Z^H
    //       ///         Z is applied to QZ
    //       {
    //         int r_val = 0;
    //         real_type *ww = w_now; w_now += (5*m); wlen_now -= (5*m);
    //         assert( (wlen_now >= 0) && "Eigendecomposition: Schur
    //         decomposition workspace ww allocation fails"); do {
    //           const bool restart = (r_val < 0);
    //           r_val = SerialSchurInternal::invoke(m,
    //                                               A, as0, as1,
    //                                               QZ, qzs0, qzs1,
    //                                               ww, 5*m,
    //                                               restart);
    //         } while (r_val < 0 && false);
    //         // for a testing purpose, we run the Schur decomposition with a
    //         finite number of Francis iterations w_now -= (5*m); wlen_now +=
    //         (5*m);
    //       }

    //       /// Step 3: Extract iigenvalues and eigenvectors from T = V S V^-1
    //       ///
    //       {
    //         /// extract eigenvalues
    //         real_type *AA = A-as1;
    //         int *blks = (int*)w_now; w_now += m; wlen_now -= m;
    //         assert( (wlen_now >= 0) && "Eigendecomposition: Eigenvector
    //         workspace blks allocation fails");

    //         {
    //           int i=0;
    //           for (;i<(m-1);) {
    //             const real_type subdiag = ats::abs(AA[(i+1)*as]);
    //             const real_type diag = A[i*as];
    //             if (subdiag < tol) {
    //               er[i*ers] = diag;
    //               ei[i*eis] = zero;
    //               blks[i] = 1;
    //               i+=1;
    //             } else {
    //               const real_type offdiag = ats::abs(A[i*as+as1]);
    //               const real_type sqrt_mult_suboffdiags =
    //               ats::sqrt(subdiag*offdiag); er[(i  )*ers] =  diag;
    //               er[(i+1)*ers] =  diag;
    //               ei[(i  )*eis] =  sqrt_mult_suboffdiags;
    //               ei[(i+1)*eis] = -sqrt_mult_suboffdiags;
    //               blks[i  ] = 2;
    //               blks[i+1] = 2; /// consider backward iteration
    //               i+=2;
    //             }
    //           }
    //           if (i<m) { /// take care the remainder
    //             er[i*ers] = A[i*as];
    //             ei[i*eis] = zero;
    //             blks[i] = 1;
    //           }
    //         }

    //         {
    //           real_type *V = w_now; w_now += (m*m); wlen_now -= (m*m);
    //           assert( (wlen_now >= 0) && "Eigendecomposition: Eigenvector
    //           workspace V allocation fails");

    //           const int vs0 = 1, vs1 = m;
    //           real_type *ww = w_now; w_now += 2*m; wlen_now -= 2*m;
    //           assert( (wlen_now >= 0) && "Eigendecomposition: Eigenvector
    //           workspace w allocation fails");

    //           /// Right eigenvectors V
    //           if (is_UR) {
    //             SerialRightEigenvectorFromSchurInternal
    //               ::invoke(m,
    //                        A, as0, as1,
    //                        V, vs0, vs1,
    //                        ww,
    //                        blks);

    //             /// QZ V
    //             SerialGemmInternal<Algo::Gemm::Unblocked>::
    //               invoke(m, m, m,
    //                      one,
    //                      QZ, qzs0, qzs1,
    //                      V, vs0, vs1,
    //                      zero,
    //                      UR, urs0, urs1);
    //             int j=0;
    //             for (;j<m;) {
    //               if (ats::abs(ei[j*eis]) < tol) {
    //                 /// a real eigenvalue
    //                 SerialNormalizeInternal::invoke(m,
    //                                                 UR+j*urs1, urs0);
    //                 j+=1;
    //               } else {
    //                 /// a complex pair of eigenvalues
    //                 SerialNormalizeInternal::invoke(m,
    //                                                 UR+(j  )*urs1, urs0,
    //                                                 UR+(j+1)*urs1, urs0);
    //                 j+=2;
    //               }
    //             }
    //           }

    //           /// Left eigenvectors V stores V
    //           if (is_UL) {
    //             SerialLeftEigenvectorFromSchurInternal
    //               ::invoke(m,
    //                        A, as0, as1,
    //                        V, vs0, vs1,
    //                        ww,
    //                        blks);

    //             /// V^-1 (QZ)^H
    //             SerialGemmInternal<Algo::Gemm::Unblocked>::
    //               invoke(m, m, m,
    //                      one,
    //                      V, vs0, vs1,
    //                      QZ, qzs1, qzs0,
    //                      zero,
    //                      UL, uls0, uls1);

    //             int i=0;
    //             for (;i<m;) {
    //               /// normalize row vectors
    //               if (ats::abs(ei[i*eis]) < tol) {
    //                 /// a real eigenvalue
    //                 SerialNormalizeInternal::invoke(m,
    //                                                 UL+i*uls0, uls1);
    //                 i+=1;
    //               } else {
    //                 /// a complex pair of eigenvalues
    //                 SerialNormalizeInternal::invoke(m,
    //                                                 UL+(i  )*uls0, uls1,
    //                                                 UL+(i+1)*uls0, uls1);
    //                 i+=2;
    //               }
    //             }
    //           }
    //           w_now -= (m*m+2*m); wlen_now += (m*m+2*m);
    //         }
    //         // deallocate blks
    //         w_now -= m; wlen_now += m;
    //       }
    //       /// deallocate QZ
    //       w_now -= (m*m); wlen_now += (m*m);

    //       assert( (w == w_now) && "Eigendecomposition: workspace tracking
    //       fails"); assert( (wlen == wlen_now) && "Eigendecomposition:
    //       workspace counting fails");
    return 0;
  }

  template <typename RealType>
  inline static int host_invoke(const int m, RealType* A, const int as0, const int as1, RealType* er, const int ers,
                                RealType* ei, const int eis, RealType* UL, const int uls0, const int uls1, RealType* UR,
                                const int urs0, const int urs1, RealType* w, const int wlen) {
#if defined(__KOKKOSBATCHED_ENABLE_LAPACKE__) || defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL__)
    int matrix_layout(0), lda(0), uls(0), urs(0);
    if (as0 == 1) {
      assert(uls0 == 1 && "UL is not column major");
      assert(urs0 == 1 && "UL is not column major");

      matrix_layout = LAPACK_COL_MAJOR;
      lda           = as1;
      uls           = uls1;
      urs           = urs1;
    }

    if (as1 == 1) {
      matrix_layout = LAPACK_ROW_MAJOR;
      assert(uls1 == 1 && "UR is not row major");
      assert(urs1 == 1 && "UR is not row major");
      lda = as0;
      uls = uls0;
      urs = urs0;
    }
    assert(matrix_layout != 0 && "Either stride of A is not unit");
    if (std::is_same<RealType, float>::value) {
      LAPACKE_sgeev(matrix_layout, 'V', 'V', m, (float*)A, lda, (float*)er, (float*)ei, (float*)UL, uls, (float*)UR,
                    urs);
    } else if (std::is_same<RealType, double>::value) {
      LAPACKE_dgeev(matrix_layout, 'V', 'V', m, (double*)A, lda, (double*)er, (double*)ei, (double*)UL, uls,
                    (double*)UR, urs);
    } else {
      // no complex is needed for this moment
      assert(false && "complex type is not supported");
    }
#else
    device_invoke(m, A, as0, as1, er, ers, ei, eis, UL, uls0, uls1, UR, urs0, urs1, w, wlen);

#endif
    return 0;
  }

  template <typename RealType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, RealType* A, const int as0, const int as1, RealType* er,
                                           const int ers, RealType* ei, const int eis, RealType* UL, const int uls0,
                                           const int uls1, RealType* UR, const int urs0, const int urs1, RealType* w,
                                           const int wlen) {
    KOKKOS_IF_ON_HOST((host_invoke(m, A, as0, as1, er, ers, ei, eis, UL, uls0, uls1, UR, urs0, urs1, w, wlen);))
    KOKKOS_IF_ON_DEVICE((device_invoke(m, A, as0, as1, er, ers, ei, eis, UL, uls0, uls1, UR, urs0, urs1, w, wlen);))
    return 0;
  }
};

}  // namespace KokkosBatched

#endif
