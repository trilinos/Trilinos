// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \file KokkosLapack_Host_tpl.cpp
/// \brief LAPACK wrapper for host tpls
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_config.h"
#include "KokkosLapack_Host_tpl.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
#include <Accelerate/Accelerate.h>
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK)

/// Fortran headers
extern "C" {

///
/// Gesv
///

void F77_BLAS_MANGLE(sgesv, SGESV)(int*, int*, float*, int*, int*, float*, int*, int*);
void F77_BLAS_MANGLE(dgesv, DGESV)(int*, int*, double*, int*, int*, double*, int*, int*);
void F77_BLAS_MANGLE(cgesv, CGESV)(int*, int*, std::complex<float>*, int*, int*, std::complex<float>*, int*, int*);
void F77_BLAS_MANGLE(zgesv, ZGESV)(int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);

///
/// Gesvd
///

void F77_BLAS_MANGLE(sgesvd, SGESVD)(const char*, const char*, const int*, const int*, float*, const int*, float*,
                                     float*, const int*, float*, const int*, float*, int*, int*);
void F77_BLAS_MANGLE(dgesvd, DGESVD)(const char*, const char*, const int*, const int*, double*, const int*, double*,
                                     double*, const int*, double*, const int*, double*, int*, int*);
void F77_BLAS_MANGLE(cgesvd, CGESVD)(const char*, const char*, const int*, const int*, std::complex<float>*, const int*,
                                     float*, std::complex<float>*, const int*, std::complex<float>*, const int*,
                                     std::complex<float>*, int*, float*, int*);
void F77_BLAS_MANGLE(zgesvd, ZGESVD)(const char*, const char*, const int*, const int*, std::complex<double>*,
                                     const int*, double*, std::complex<double>*, const int*, std::complex<double>*,
                                     const int*, std::complex<double>*, int*, double*, int*);

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
void F77_BLAS_MANGLE(strtri, STRTRI)(const char*, const char*, int*, const float*, int*, int*);
void F77_BLAS_MANGLE(dtrtri, DTRTRI)(const char*, const char*, int*, const double*, int*, int*);
void F77_BLAS_MANGLE(ctrtri, CTRTRI)(const char*, const char*, int*, const std::complex<float>*, int*, int*);
void F77_BLAS_MANGLE(ztrtri, ZTRTRI)(const char*, const char*, int*, const std::complex<double>*, int*, int*);

///
/// Geqrf
///

void F77_BLAS_MANGLE(sgeqrf, SGEQRF)(const int*, const int*, float*, const int*, float*, float*, int*, int*);
void F77_BLAS_MANGLE(dgeqrf, DGEQRF)(const int*, const int*, double*, const int*, double*, double*, int*, int*);
void F77_BLAS_MANGLE(cgeqrf, CGEQRF)(const int*, const int*, std::complex<float>*, const int*, std::complex<float>*,
                                     std::complex<float>*, int*, int*);
void F77_BLAS_MANGLE(zgeqrf, ZGEQRF)(const int*, const int*, std::complex<double>*, const int*, std::complex<double>*,
                                     std::complex<double>*, int*, int*);

///
/// {Un,Or}mqr
///
void F77_BLAS_MANGLE(sormqr, SORMQR)(const char*, const char*, const int*, const int*, const int*, float*, const int*,
                                     float*, float*, const int*, float*, int*, int*);
void F77_BLAS_MANGLE(dormqr, DORMQR)(const char*, const char*, const int*, const int*, const int*, double*, const int*,
                                     double*, double*, const int*, double*, int*, int*);
void F77_BLAS_MANGLE(cunmqr, CUNMQR)(const char*, const char*, const int*, const int*, const int*, std::complex<float>*,
                                     const int*, std::complex<float>*, std::complex<float>*, const int*,
                                     std::complex<float>*, int*, int*);
void F77_BLAS_MANGLE(zunmqr, ZUNMQR)(const char*, const char*, const int*, const int*, const int*,
                                     std::complex<double>*, const int*, std::complex<double>*, std::complex<double>*,
                                     const int*, std::complex<double>*, int*, int*);

///
/// {Un,Or}gqr
///
void F77_BLAS_MANGLE(sorgqr, SORGQR)(const int*, const int*, const int*, float*, const int*, float*, float*, int*,
                                     int*);
void F77_BLAS_MANGLE(dorgqr, DORGQR)(const int*, const int*, const int*, double*, const int*, double*, double*, int*,
                                     int*);
void F77_BLAS_MANGLE(cungqr, CUNGQR)(const int*, const int*, const int*, std::complex<float>*, const int*,
                                     std::complex<float>*, std::complex<float>*, int*, int*);
void F77_BLAS_MANGLE(zungqr, ZUNGQR)(const int*, const int*, const int*, std::complex<double>*, const int*,
                                     std::complex<double>*, std::complex<double>*, int*, int*);

///
/// Potrf
///

void F77_BLAS_MANGLE(spotrf, SPOTRF)(const char*, const int*, float*, const int*, int*);
void F77_BLAS_MANGLE(dpotrf, DPOTRF)(const char*, const int*, double*, const int*, int*);
void F77_BLAS_MANGLE(cpotrf, CPOTRF)(const char*, const int*, std::complex<float>*, const int*, int*);
void F77_BLAS_MANGLE(zpotrf, ZPOTRF)(const char*, const int*, std::complex<double>*, const int*, int*);

///
/// Potrs
///

void F77_BLAS_MANGLE(spotrs, SPOTRS)(const char*, const int*, const int*, const float*, const int*, float*, const int*,
                                     int*);
void F77_BLAS_MANGLE(dpotrs, DPOTRS)(const char*, const int*, const int*, const double*, const int*, double*,
                                     const int*, int*);
void F77_BLAS_MANGLE(cpotrs, CPOTRS)(const char*, const int*, const int*, const std::complex<float>*, const int*,
                                     std::complex<float>*, const int*, int*);
void F77_BLAS_MANGLE(zpotrs, ZPOTRS)(const char*, const int*, const int*, const std::complex<double>*, const int*,
                                     std::complex<double>*, const int*, int*);
}

#define F77_FUNC_SGESV F77_BLAS_MANGLE(sgesv, SGESV)
#define F77_FUNC_DGESV F77_BLAS_MANGLE(dgesv, DGESV)
#define F77_FUNC_CGESV F77_BLAS_MANGLE(cgesv, CGESV)
#define F77_FUNC_ZGESV F77_BLAS_MANGLE(zgesv, ZGESV)

#define F77_FUNC_SGESVD F77_BLAS_MANGLE(sgesvd, SGESVD)
#define F77_FUNC_DGESVD F77_BLAS_MANGLE(dgesvd, DGESVD)
#define F77_FUNC_CGESVD F77_BLAS_MANGLE(cgesvd, CGESVD)
#define F77_FUNC_ZGESVD F77_BLAS_MANGLE(zgesvd, ZGESVD)

#define F77_FUNC_STRTRI F77_BLAS_MANGLE(strtri, STRTRI)
#define F77_FUNC_DTRTRI F77_BLAS_MANGLE(dtrtri, DTRTRI)
#define F77_FUNC_CTRTRI F77_BLAS_MANGLE(ctrtri, CTRTRI)
#define F77_FUNC_ZTRTRI F77_BLAS_MANGLE(ztrtri, ZTRTRI)

#define F77_FUNC_SGEQRF F77_BLAS_MANGLE(sgeqrf, SGEQRF)
#define F77_FUNC_DGEQRF F77_BLAS_MANGLE(dgeqrf, DGEQRF)
#define F77_FUNC_CGEQRF F77_BLAS_MANGLE(cgeqrf, CGEQRF)
#define F77_FUNC_ZGEQRF F77_BLAS_MANGLE(zgeqrf, ZGEQRF)

#define F77_FUNC_SORMQR F77_BLAS_MANGLE(sormqr, SORMQR)
#define F77_FUNC_DORMQR F77_BLAS_MANGLE(dormqr, DORMQR)
#define F77_FUNC_CUNMQR F77_BLAS_MANGLE(cunmqr, CUNMQR)
#define F77_FUNC_ZUNMQR F77_BLAS_MANGLE(zunmqr, ZUNMQR)

#define F77_FUNC_SORGQR F77_BLAS_MANGLE(sorgqr, SORGQR)
#define F77_FUNC_DORGQR F77_BLAS_MANGLE(dorgqr, DORGQR)
#define F77_FUNC_CUNGQR F77_BLAS_MANGLE(cungqr, CUNGQR)
#define F77_FUNC_ZUNGQR F77_BLAS_MANGLE(zungqr, ZUNGQR)

#define F77_FUNC_SPOTRF F77_BLAS_MANGLE(spotrf, SPOTRF)
#define F77_FUNC_DPOTRF F77_BLAS_MANGLE(dpotrf, DPOTRF)
#define F77_FUNC_CPOTRF F77_BLAS_MANGLE(cpotrf, CPOTRF)
#define F77_FUNC_ZPOTRF F77_BLAS_MANGLE(zpotrf, ZPOTRF)

#define F77_FUNC_SPOTRS F77_BLAS_MANGLE(spotrs, SPOTRS)
#define F77_FUNC_DPOTRS F77_BLAS_MANGLE(dpotrs, DPOTRS)
#define F77_FUNC_CPOTRS F77_BLAS_MANGLE(cpotrs, CPOTRS)
#define F77_FUNC_ZPOTRS F77_BLAS_MANGLE(zpotrs, ZPOTRS)

#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
namespace KokkosLapack {
namespace Impl {

///
/// float
///

template <>
void HostLapack<float>::gesv(int n, int rhs, float* a, int lda, int* ipiv, float* b, int ldb, int info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  sgesv_(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
#else
  F77_FUNC_SGESV(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
#endif
}

template <>
void HostLapack<float>::gesvd(const char jobu, const char jobvt, const int m, const int n, float* a, const int lda,
                              float* s, float* u, const int ldu, float* vt, const int ldvt, float* work, int lwork,
                              float* /*rwork*/, int info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  sgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
#else
  F77_FUNC_SGESVD(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
#endif
}

template <>
int HostLapack<float>::trtri(const char uplo, const char diag, int n, const float* a, int lda) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  strtri_(&uplo, &diag, &n, const_cast<float*>(a), &lda, &info);
#else
  F77_FUNC_STRTRI(&uplo, &diag, &n, a, &lda, &info);
#endif
  return info;
}

template <>
void HostLapack<float>::geqrf(const int m, const int n, float* a, const int lda, float* tau, float* work, int lwork,
                              int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  sgeqrf_(&m, &n, a, &lda, tau, work, &lwork, info);
#else
  F77_FUNC_SGEQRF(&m, &n, a, &lda, tau, work, &lwork, info);
#endif
}
template <>
void HostLapack<float>::gemqr(const char side, const char trans, const int m, const int n, const int k, float* a,
                              const int lda, float* tau, float* c, const int ldc, float* work, int lwork, int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  sormqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info);
#else
  F77_FUNC_SORMQR(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info);
#endif
}
template <>
void HostLapack<float>::gegqr(const int m, const int n, const int k, float* a, const int lda, float* tau, float* work,
                              int lwork, int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  sorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
#else
  F77_FUNC_SORGQR(&m, &n, &k, a, &lda, tau, work, &lwork, info);
#endif
}

template <>
int HostLapack<float>::potrf(const char uplo, const int n, float* a, const int lda) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  spotrf_(&uplo, &n, a, &lda, &info);
#else
  F77_FUNC_SPOTRF(&uplo, &n, a, &lda, &info);
#endif
  return info;
}
template <>
int HostLapack<float>::potrs(const char uplo, const int n, const int nrhs, const float* a, const int lda, float* b,
                             const int ldb) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  spotrs_(&uplo, &n, &nrhs, const_cast<float*>(a), &lda, b, &ldb, &info);
#else
  F77_FUNC_SPOTRS(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
#endif
  return info;
}

///
/// double
///

template <>
void HostLapack<double>::gesv(int n, int rhs, double* a, int lda, int* ipiv, double* b, int ldb, int info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  dgesv_(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
#else
  F77_FUNC_DGESV(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
#endif
}

template <>
void HostLapack<double>::gesvd(const char jobu, const char jobvt, const int m, const int n, double* a, const int lda,
                               double* s, double* u, const int ldu, double* vt, const int ldvt, double* work, int lwork,
                               double* /*rwork*/, int info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
#else
  F77_FUNC_DGESVD(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
#endif
}

template <>
int HostLapack<double>::trtri(const char uplo, const char diag, int n, const double* a, int lda) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  dtrtri_(&uplo, &diag, &n, const_cast<double*>(a), &lda, &info);
#else
  F77_FUNC_DTRTRI(&uplo, &diag, &n, a, &lda, &info);
#endif
  return info;
}

template <>
void HostLapack<double>::geqrf(const int m, const int n, double* a, const int lda, double* tau, double* work, int lwork,
                               int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, info);
#else
  F77_FUNC_DGEQRF(&m, &n, a, &lda, tau, work, &lwork, info);
#endif
}

template <>
int HostLapack<double>::potrf(const char uplo, const int n, double* a, const int lda) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  dpotrf_(&uplo, &n, a, &lda, &info);
#else
  F77_FUNC_DPOTRF(&uplo, &n, a, &lda, &info);
#endif
  return info;
}

template <>
int HostLapack<double>::potrs(const char uplo, const int n, const int nrhs, const double* a, const int lda, double* b,
                              const int ldb) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  dpotrs_(&uplo, &n, &nrhs, const_cast<double*>(a), &lda, b, &ldb, &info);
#else
  F77_FUNC_DPOTRS(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
#endif
  return info;
}

template <>
void HostLapack<double>::gemqr(const char side, const char trans, const int m, const int n, const int k, double* a,
                               const int lda, double* tau, double* c, const int ldc, double* work, int lwork,
                               int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  dormqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info);
#else
  F77_FUNC_DORMQR(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info);
#endif
}
template <>
void HostLapack<double>::gegqr(const int m, const int n, const int k, double* a, const int lda, double* tau,
                               double* work, int lwork, int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
#else
  F77_FUNC_DORGQR(&m, &n, &k, a, &lda, tau, work, &lwork, info);
#endif
}

///
/// std::complex<float>
///

template <>
void HostLapack<std::complex<float>>::gesv(int n, int rhs, std::complex<float>* a, int lda, int* ipiv,
                                           std::complex<float>* b, int ldb, int info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cgesv_(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
#else
  F77_FUNC_CGESV(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
#endif
}
template <>
void HostLapack<std::complex<float>>::gesvd(const char jobu, const char jobvt, const int m, const int n,
                                            std::complex<float>* a, const int lda, float* s, std::complex<float>* u,
                                            const int ldu, std::complex<float>* vt, const int ldvt,
                                            std::complex<float>* work, int lwork, float* rwork, int info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);
#else
  F77_FUNC_CGESVD(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);
#endif
}
template <>
int HostLapack<std::complex<float>>::trtri(const char uplo, const char diag, int n, const std::complex<float>* a,
                                           int lda) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  ctrtri_(&uplo, &diag, &n, const_cast<std::complex<float>*>(a), &lda, &info);
#else
  F77_FUNC_CTRTRI(&uplo, &diag, &n, a, &lda, &info);
#endif
  return info;
}

template <>
void HostLapack<std::complex<float>>::geqrf(const int m, const int n, std::complex<float>* a, const int lda,
                                            std::complex<float>* tau, std::complex<float>* work, int lwork, int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cgeqrf_(&m, &n, a, &lda, tau, work, &lwork, info);
#else
  F77_FUNC_CGEQRF(&m, &n, a, &lda, tau, work, &lwork, info);
#endif
}
template <>
void HostLapack<std::complex<float>>::gemqr(const char side, const char trans, const int m, const int n, const int k,
                                            std::complex<float>* a, const int lda, std::complex<float>* tau,
                                            std::complex<float>* c, const int ldc, std::complex<float>* work, int lwork,
                                            int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cunmqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info);
#else
  F77_FUNC_CUNMQR(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info);
#endif
}
template <>
void HostLapack<std::complex<float>>::gegqr(const int m, const int n, const int k, std::complex<float>* a,
                                            const int lda, std::complex<float>* tau, std::complex<float>* work,
                                            int lwork, int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cungqr_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
#else
  F77_FUNC_CUNGQR(&m, &n, &k, a, &lda, tau, work, &lwork, info);
#endif
}

template <>
int HostLapack<std::complex<float>>::potrf(const char uplo, const int n, std::complex<float>* a, const int lda) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cpotrf_(&uplo, &n, a, &lda, &info);
#else
  F77_FUNC_CPOTRF(&uplo, &n, a, &lda, &info);
#endif
  return info;
}
template <>
int HostLapack<std::complex<float>>::potrs(const char uplo, const int n, const int nrhs, const std::complex<float>* a,
                                           const int lda, std::complex<float>* b, const int ldb) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cpotrs_(&uplo, &n, &nrhs, const_cast<std::complex<float>*>(a), &lda, b, &ldb, &info);
#else
  F77_FUNC_CPOTRS(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
#endif
  return info;
}

///
/// std::complex<double>
///

template <>
void HostLapack<std::complex<double>>::gesv(int n, int rhs, std::complex<double>* a, int lda, int* ipiv,
                                            std::complex<double>* b, int ldb, int info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  zgesv_(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
#else
  F77_FUNC_ZGESV(&n, &rhs, a, &lda, ipiv, b, &ldb, &info);
#endif
}
template <>
void HostLapack<std::complex<double>>::gesvd(const char jobu, const char jobvt, const int m, const int n,
                                             std::complex<double>* a, const int lda, double* s, std::complex<double>* u,
                                             const int ldu, std::complex<double>* vt, const int ldvt,
                                             std::complex<double>* work, int lwork, double* rwork, int info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  zgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);
#else
  F77_FUNC_ZGESVD(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);
#endif
}
template <>
int HostLapack<std::complex<double>>::trtri(const char uplo, const char diag, int n, const std::complex<double>* a,
                                            int lda) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  ztrtri_(&uplo, &diag, &n, const_cast<std::complex<double>*>(a), &lda, &info);
#else
  F77_FUNC_ZTRTRI(&uplo, &diag, &n, a, &lda, &info);
#endif
  return info;
}

template <>
void HostLapack<std::complex<double>>::geqrf(const int m, const int n, std::complex<double>* a, const int lda,
                                             std::complex<double>* tau, std::complex<double>* work, int lwork,
                                             int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  zgeqrf_(&m, &n, a, &lda, tau, work, &lwork, info);
#else
  F77_FUNC_ZGEQRF(&m, &n, a, &lda, tau, work, &lwork, info);
#endif
}
template <>
void HostLapack<std::complex<double>>::gemqr(const char side, const char trans, const int m, const int n, const int k,
                                             std::complex<double>* a, const int lda, std::complex<double>* tau,
                                             std::complex<double>* c, const int ldc, std::complex<double>* work,
                                             int lwork, int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  zunmqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info);
#else
  F77_FUNC_ZUNMQR(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info);
#endif
}
template <>
void HostLapack<std::complex<double>>::gegqr(const int m, const int n, const int k, std::complex<double>* a,
                                             const int lda, std::complex<double>* tau, std::complex<double>* work,
                                             int lwork, int* info) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  zungqr_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
#else
  F77_FUNC_ZUNGQR(&m, &n, &k, a, &lda, tau, work, &lwork, info);
#endif
}

template <>
int HostLapack<std::complex<double>>::potrf(const char uplo, const int n, std::complex<double>* a, const int lda) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  zpotrf_(&uplo, &n, a, &lda, &info);
#else
  F77_FUNC_ZPOTRF(&uplo, &n, a, &lda, &info);
#endif
  return info;
}
template <>
int HostLapack<std::complex<double>>::potrs(const char uplo, const int n, const int nrhs, const std::complex<double>* a,
                                            const int lda, std::complex<double>* b, const int ldb) {
  int info = 0;
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  zpotrs_(&uplo, &n, &nrhs, const_cast<std::complex<double>*>(a), &lda, b, &ldb, &info);
#else
  F77_FUNC_ZPOTRS(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
#endif
  return info;
}

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK OR KOKKOSKERNELS_ENABLE_TPL_ACCELERATE
