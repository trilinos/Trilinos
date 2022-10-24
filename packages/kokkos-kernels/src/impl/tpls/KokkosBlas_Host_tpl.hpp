/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSBLAS_HOST_TPL_HPP_
#define KOKKOSBLAS_HOST_TPL_HPP_

/// \file  KokkosBlas_Host_tpl.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_config.h"
#include "Kokkos_ArithTraits.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)

namespace KokkosBlas {
namespace Impl {

template <typename T>
struct HostBlas {
  typedef Kokkos::ArithTraits<T> ats;
  typedef typename ats::mag_type mag_type;

  static void scal(int n, const T alpha,
                   /* */ T *x, int x_inc);

  static int iamax(int n, const T *x, int x_inc);

  static mag_type nrm2(int n, const T *x, int x_inc);

  static mag_type asum(int n, const T *x, int x_inc);

  static T dot(int n, const T *x, int x_inc, const T *y, int y_inc);

  static void axpy(int n, const T alpha, const T *x, int x_inc,
                   /* */ T *y, int y_inc);

  static void gemv(const char trans, int m, int n, const T alpha, const T *a,
                   int lda, const T *b, int ldb, const T beta,
                   /* */ T *c, int ldc);

  static void trsv(const char uplo, const char transa, const char diag, int m,
                   const T *a, int lda,
                   /* */ T *b, int ldb);

  static void gemm(const char transa, const char transb, int m, int n, int k,
                   const T alpha, const T *a, int lda, const T *b, int ldb,
                   const T beta,
                   /* */ T *c, int ldc);

  static void herk(const char transa, const char transb, int n, int k,
                   const T alpha, const T *a, int lda, const T beta,
                   /* */ T *c, int ldc);

  static void trmm(const char side, const char uplo, const char transa,
                   const char diag, int m, int n, const T alpha, const T *a,
                   int lda,
                   /* */ T *b, int ldb);

  static void trsm(const char side, const char uplo, const char transa,
                   const char diag, int m, int n, const T alpha, const T *a,
                   int lda,
                   /* */ T *b, int ldb);

  static void gesv(int n, int rhs, T *a, int lda, int *ipiv, T *b, int ldb,
                   int info);

  static int trtri(const char uplo, const char diag, int n, const T *a,
                   int lda);
};
}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

#endif  // KOKKOSBLAS_HOST_TPL_HPP_
