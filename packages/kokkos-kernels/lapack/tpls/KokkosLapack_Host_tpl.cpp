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
/// \file KokkosLapack_Host_tpl.cpp
/// \brief LAPACK wrapper for host tpls
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_config.h"
#include "KokkosLapack_Host_tpl.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK)

/// Fortran headers
extern "C" {

///
/// Gesv
///

void F77_BLAS_MANGLE(sgesv, SGESV)(int*, int*, float*, int*, int*, float*, int*,
                                   int*);
void F77_BLAS_MANGLE(dgesv, DGESV)(int*, int*, double*, int*, int*, double*,
                                   int*, int*);
void F77_BLAS_MANGLE(cgesv, CGESV)(int*, int*, std::complex<float>*, int*, int*,
                                   std::complex<float>*, int*, int*);
void F77_BLAS_MANGLE(zgesv, ZGESV)(int*, int*, std::complex<double>*, int*,
                                   int*, std::complex<double>*, int*, int*);

///
/// Trtri
///
/*
    HostLapack<float>::trtri(const char uplo, const char diag,
                           int n, const float *a, int lda) {
      int info = 0;
      F77_FUNC_STRTRI(&uplo,
                      &diag, &n,
                      a, &lda, &info);
*/
void F77_BLAS_MANGLE(strtri, STRTRI)(const char*, const char*, int*,
                                     const float*, int*, int*);
void F77_BLAS_MANGLE(dtrtri, DTRTRI)(const char*, const char*, int*,
                                     const double*, int*, int*);
void F77_BLAS_MANGLE(ctrtri, CTRTRI)(const char*, const char*, int*,
                                     const std::complex<float>*, int*, int*);
void F77_BLAS_MANGLE(ztrtri, ZTRTRI)(const char*, const char*, int*,
                                     const std::complex<double>*, int*, int*);
}

#define F77_FUNC_SGESV F77_BLAS_MANGLE(sgesv, SGESV)
#define F77_FUNC_DGESV F77_BLAS_MANGLE(dgesv, DGESV)
#define F77_FUNC_CGESV F77_BLAS_MANGLE(cgesv, CGESV)
#define F77_FUNC_ZGESV F77_BLAS_MANGLE(zgesv, ZGESV)

#define F77_FUNC_STRTRI F77_BLAS_MANGLE(strtri, STRTRI)
#define F77_FUNC_DTRTRI F77_BLAS_MANGLE(dtrtri, DTRTRI)
#define F77_FUNC_CTRTRI F77_BLAS_MANGLE(ctrtri, CTRTRI)
#define F77_FUNC_ZTRTRI F77_BLAS_MANGLE(ztrtri, ZTRTRI)

namespace KokkosLapack {
namespace Impl {

///
/// float
///

template <>
void HostLapack<float>::gesv(int n, int rhs, float* a, int lda, int* ipiv,
                             float* b, int ldb, int info) {
  F77_FUNC_SGESV(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
}
template <>
int HostLapack<float>::trtri(const char uplo, const char diag, int n,
                             const float* a, int lda) {
  int info = 0;
  F77_FUNC_STRTRI(&uplo, &diag, &n, a, &lda, &info);
  return info;
}

///
/// double
///

template <>
void HostLapack<double>::gesv(int n, int rhs, double* a, int lda, int* ipiv,
                              double* b, int ldb, int info) {
  F77_FUNC_DGESV(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
}
template <>
int HostLapack<double>::trtri(const char uplo, const char diag, int n,
                              const double* a, int lda) {
  int info = 0;
  F77_FUNC_DTRTRI(&uplo, &diag, &n, a, &lda, &info);
  return info;
}

///
/// std::complex<float>
///

template <>
void HostLapack<std::complex<float> >::gesv(int n, int rhs,
                                            std::complex<float>* a, int lda,
                                            int* ipiv, std::complex<float>* b,
                                            int ldb, int info) {
  F77_FUNC_CGESV(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
}
template <>
int HostLapack<std::complex<float> >::trtri(const char uplo, const char diag,
                                            int n, const std::complex<float>* a,
                                            int lda) {
  int info = 0;
  F77_FUNC_CTRTRI(&uplo, &diag, &n, a, &lda, &info);
  return info;
}

///
/// std::complex<double>
///

template <>
void HostLapack<std::complex<double> >::gesv(int n, int rhs,
                                             std::complex<double>* a, int lda,
                                             int* ipiv, std::complex<double>* b,
                                             int ldb, int info) {
  F77_FUNC_ZGESV(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
}
template <>
int HostLapack<std::complex<double> >::trtri(const char uplo, const char diag,
                                             int n,
                                             const std::complex<double>* a,
                                             int lda) {
  int info = 0;
  F77_FUNC_ZTRTRI(&uplo, &diag, &n, a, &lda, &info);
  return info;
}

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK
