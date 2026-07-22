// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \file KokkosBlas_Host_tpl.cpp
/// \brief BLAS wrapper for host tpls
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_config.h"
#include "KokkosBlas_Host_tpl.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)

using KokkosBlas::Impl::KK_INT;

#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
#include <Accelerate/Accelerate.h>

auto kk_to_accelerate_trans(const char trans) {
  if (trans == 'N' || trans == 'n') {
    return CblasNoTrans;
  } else if (trans == 'T' || trans == 't') {
    return CblasTrans;
  } else if (trans == 'C' || trans == 'c') {
    return CblasConjTrans;
  } else {
    throw std::invalid_argument("KokkosBlas[cblas]: invalid trans parameter");
  }
}

CBLAS_UPLO kk_to_accelerate_uplo(const char uplo) {
  if (uplo == 'U' || uplo == 'u') {
    return CblasUpper;
  } else if (uplo == 'L' || uplo == 'l') {
    return CblasLower;
  } else {
    throw std::invalid_argument("KokkosBlas[cblas]: invalid uplo parameter");
  }
}

CBLAS_DIAG kk_to_accelerate_diag(const char diag) {
  if (diag == 'U' || diag == 'u') {
    return CblasUnit;
  } else if (diag == 'N' || diag == 'n') {
    return CblasNonUnit;
  } else {
    throw std::invalid_argument("KokkosBlas[cblas]: invalid diag parameter");
  }
}

CBLAS_SIDE kk_to_accelerate_side(const char side) {
  if (side == 'L' || side == 'l') {
    return CblasLeft;
  } else if (side == 'R' || side == 'r') {
    return CblasRight;
  } else {
    throw std::invalid_argument("KokkosBlas[cblas]: invalid side parameter");
  }
}

#else
/// Fortran headers
extern "C" {

///
/// scal
///
void F77_BLAS_MANGLE(sscal, SSCAL)(const KK_INT* N, const float* alpha,
                                   /* */ float* x, const KK_INT* x_inc);
void F77_BLAS_MANGLE(dscal, DSCAL)(const KK_INT* N, const double* alpha,
                                   /* */ double* x, const KK_INT* x_inc);
void F77_BLAS_MANGLE(cscal, CSCAL)(const KK_INT* N, const std::complex<float>* alpha,
                                   /* */ std::complex<float>* x, const KK_INT* x_inc);
void F77_BLAS_MANGLE(zscal, ZSCAL)(const KK_INT* N, const std::complex<double>* alpha,
                                   /* */ std::complex<double>* x, const KK_INT* x_inc);

///
/// max
///
KK_INT F77_BLAS_MANGLE(isamax, ISAMAX)(const KK_INT* N, const float* x, const KK_INT* x_inc);
KK_INT F77_BLAS_MANGLE(idamax, IDAMAX)(const KK_INT* N, const double* x, const KK_INT* x_inc);
KK_INT F77_BLAS_MANGLE(icamax, ICAMAX)(const KK_INT* N, const std::complex<float>* x, const KK_INT* x_inc);
KK_INT F77_BLAS_MANGLE(izamax, IZAMAX)(const KK_INT* N, const std::complex<double>* x, const KK_INT* x_inc);

///
/// nrm2
///
float F77_BLAS_MANGLE(snrm2, SNRM2)(const KK_INT* N, const float* x, const KK_INT* x_inc);
double F77_BLAS_MANGLE(dnrm2, DNRM2)(const KK_INT* N, const double* x, const KK_INT* x_inc);
float F77_BLAS_MANGLE(scnrm2, SCNRM2)(const KK_INT* N, const std::complex<float>* x, const KK_INT* x_inc);
double F77_BLAS_MANGLE(dznrm2, DZNRM2)(const KK_INT* N, const std::complex<double>* x, const KK_INT* x_inc);

///
/// sum
///
float F77_BLAS_MANGLE(sasum, SASUM)(const KK_INT* N, const float* x, const KK_INT* x_inc);
double F77_BLAS_MANGLE(dasum, DASUM)(const KK_INT* N, const double* x, const KK_INT* x_inc);
float F77_BLAS_MANGLE(scasum, SCASUM)(const KK_INT* N, const std::complex<float>* x, const KK_INT* x_inc);
double F77_BLAS_MANGLE(dzasum, DZASUM)(const KK_INT* N, const std::complex<double>* x, const KK_INT* x_inc);

///
/// dot
///
float F77_BLAS_MANGLE(sdot, SDOT)(const KK_INT* N, const float* x, const KK_INT* x_inc, const float* y,
                                  const KK_INT* y_inc);
double F77_BLAS_MANGLE(ddot, DDOT)(const KK_INT* N, const double* x, const KK_INT* x_inc, const double* y,
                                   const KK_INT* y_inc);
#if defined(KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX)
// clang-format off
// For the return type, don't use std::complex, otherwise compiler will complain
// error: 'cdotu_' has C-linkage specified, but returns user-defined type 'std::complex' which is incompatible with C [-Werror,-Wreturn-type-c-linkage]"
// But with float _Complex, I got error: '_Complex' is a C99 extension [-Werror,-Wc99-extensions].
// So I just use a C struct.
// clang-format on
typedef struct {
  float vals[2];
} _kk_float2;
typedef struct {
  double vals[2];
} _kk_double2;

_kk_float2 F77_BLAS_MANGLE(cdotu, CDOTU)(const KK_INT* N, const std::complex<float>* x, const KK_INT* x_inc,
                                         const std::complex<float>* y, const KK_INT* y_inc);
_kk_double2 F77_BLAS_MANGLE(zdotu, ZDOTU)(const KK_INT* N, const std::complex<double>* x, const KK_INT* x_inc,
                                          const std::complex<double>* y, const KK_INT* y_inc);
_kk_float2 F77_BLAS_MANGLE(cdotc, CDOTC)(const KK_INT* N, const std::complex<float>* x, const KK_INT* x_inc,
                                         const std::complex<float>* y, const KK_INT* y_inc);
_kk_double2 F77_BLAS_MANGLE(zdotc, ZDOTC)(const KK_INT* N, const std::complex<double>* x, const KK_INT* x_inc,
                                          const std::complex<double>* y, const KK_INT* y_inc);
#else
void F77_BLAS_MANGLE(cdotu, CDOTU)(std::complex<float>* res, const KK_INT* N, const std::complex<float>* x,
                                   const KK_INT* x_inc, const std::complex<float>* y, const KK_INT* y_inc);
void F77_BLAS_MANGLE(zdotu, ZDOTU)(std::complex<double>* res, const KK_INT* N, const std::complex<double>* x,
                                   const KK_INT* x_inc, const std::complex<double>* y, const KK_INT* y_inc);
void F77_BLAS_MANGLE(cdotc, CDOTC)(std::complex<float>* res, const KK_INT* N, const std::complex<float>* x,
                                   const KK_INT* x_inc, const std::complex<float>* y, const KK_INT* y_inc);
void F77_BLAS_MANGLE(zdotc, ZDOTC)(std::complex<double>* res, const KK_INT* N, const std::complex<double>* x,
                                   const KK_INT* x_inc, const std::complex<double>* y, const KK_INT* y_inc);
#endif

///
/// axpy
///
void F77_BLAS_MANGLE(saxpy, SAXPY)(const KK_INT* N, const float* alpha, const float* x, const KK_INT* x_inc,
                                   /* */ float* y, const KK_INT* y_inc);
void F77_BLAS_MANGLE(daxpy, DAXPY)(const KK_INT* N, const double* alpha, const double* x, const KK_INT* x_inc,
                                   /* */ double* y, const KK_INT* y_inc);
void F77_BLAS_MANGLE(caxpy, CAXPY)(const KK_INT* N, const std::complex<float>* alpha, const std::complex<float>* x,
                                   const KK_INT* x_inc,
                                   /* */ std::complex<float>* y, const KK_INT* y_inc);
void F77_BLAS_MANGLE(zaxpy, ZAXPY)(const KK_INT* N, const std::complex<double>* alpha, const std::complex<double>* x,
                                   const KK_INT* x_inc,
                                   /* */ std::complex<double>* y, const KK_INT* y_inc);

///
/// rot
///
void F77_BLAS_MANGLE(srot, SROT)(KK_INT const* N, float* X, KK_INT const* incx, float* Y, KK_INT const* incy, float* c,
                                 float* s);
void F77_BLAS_MANGLE(drot, DROT)(KK_INT const* N, double* X, KK_INT const* incx, double* Y, KK_INT const* incy,
                                 double* c, double* s);
void F77_BLAS_MANGLE(crot, CROT)(KK_INT const* N, std::complex<float>* X, KK_INT const* incx, std::complex<float>* Y,
                                 KK_INT const* incy, float* c, std::complex<float>* s);
void F77_BLAS_MANGLE(zrot, ZROT)(KK_INT const* N, std::complex<double>* X, KK_INT const* incx, std::complex<double>* Y,
                                 KK_INT const* incy, double* c, std::complex<double>* s);

///
/// rotg
///
void F77_BLAS_MANGLE(srotg, SROTG)(float* a, float* b, float* c, float* s);
void F77_BLAS_MANGLE(drotg, DROTG)(double* a, double* b, double* c, double* s);
void F77_BLAS_MANGLE(crotg, CROTG)(std::complex<float>* a, std::complex<float>* b, float* c, std::complex<float>* s);
void F77_BLAS_MANGLE(zrotg, ZROTG)(std::complex<double>* a, std::complex<double>* b, double* c,
                                   std::complex<double>* s);

///
/// rotm
///
void F77_BLAS_MANGLE(srotm, SROTM)(const KK_INT* n, float* X, const KK_INT* incx, float* Y, const KK_INT* incy,
                                   float const* param);
void F77_BLAS_MANGLE(drotm, DROTM)(const KK_INT* n, double* X, const KK_INT* incx, double* Y, const KK_INT* incy,
                                   double const* param);

///
/// rotmg
///
void F77_BLAS_MANGLE(srotmg, SROTMG)(float* d1, float* d2, float* x1, const float* y1, float* param);
void F77_BLAS_MANGLE(drotmg, DROTMG)(double* d1, double* d2, double* x1, const double* y1, double* param);

///
/// swap
///
void F77_BLAS_MANGLE(sswap, SSWAP)(KK_INT const* N, float* X, KK_INT const* incx, float* Y, KK_INT const* incy);
void F77_BLAS_MANGLE(dswap, DSWAP)(KK_INT const* N, double* X, KK_INT const* incx, double* Y, KK_INT const* incy);
void F77_BLAS_MANGLE(cswap, CSWAP)(KK_INT const* N, std::complex<float>* X, KK_INT const* incx, std::complex<float>* Y,
                                   KK_INT const* incy);
void F77_BLAS_MANGLE(zswap, ZSWAP)(KK_INT const* N, std::complex<double>* X, KK_INT const* incx,
                                   std::complex<double>* Y, KK_INT const* incy);

///
/// Gemv
///
void F77_BLAS_MANGLE(sgemv, SGEMV)(const char*, KK_INT*, KK_INT*, const float*, const float*, KK_INT*, const float*,
                                   KK_INT*, const float*,
                                   /* */ float*, KK_INT*);
void F77_BLAS_MANGLE(dgemv, DGEMV)(const char*, KK_INT*, KK_INT*, const double*, const double*, KK_INT*, const double*,
                                   KK_INT*, const double*,
                                   /* */ double*, KK_INT*);
void F77_BLAS_MANGLE(cgemv, CGEMV)(const char*, KK_INT*, KK_INT*, const std::complex<float>*,
                                   const std::complex<float>*, KK_INT*, const std::complex<float>*, KK_INT*,
                                   const std::complex<float>*,
                                   /* */ std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(zgemv, ZGEMV)(const char*, KK_INT*, KK_INT*, const std::complex<double>*,
                                   const std::complex<double>*, KK_INT*, const std::complex<double>*, KK_INT*,
                                   const std::complex<double>*,
                                   /* */ std::complex<double>*, KK_INT*);

///
/// Ger
///
void F77_BLAS_MANGLE(sger, SGER)(KK_INT*, KK_INT*, const float*, const float*, KK_INT*, const float*, KK_INT*, float*,
                                 KK_INT*);
void F77_BLAS_MANGLE(dger, DGER)(KK_INT*, KK_INT*, const double*, const double*, KK_INT*, const double*, KK_INT*,
                                 double*, KK_INT*);
void F77_BLAS_MANGLE(cgeru, CGERU)(KK_INT*, KK_INT*, const std::complex<float>*, const std::complex<float>*, KK_INT*,
                                   const std::complex<float>*, KK_INT*, std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(zgeru, ZGERU)(KK_INT*, KK_INT*, const std::complex<double>*, const std::complex<double>*, KK_INT*,
                                   const std::complex<double>*, KK_INT*, std::complex<double>*, KK_INT*);
void F77_BLAS_MANGLE(cgerc, CGERC)(KK_INT*, KK_INT*, const std::complex<float>*, const std::complex<float>*, KK_INT*,
                                   const std::complex<float>*, KK_INT*, std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(zgerc, ZGERC)(KK_INT*, KK_INT*, const std::complex<double>*, const std::complex<double>*, KK_INT*,
                                   const std::complex<double>*, KK_INT*, std::complex<double>*, KK_INT*);

///
/// Syr
///
void F77_BLAS_MANGLE(ssyr, SSYR)(const char*, KK_INT*, const float*, const float*, KK_INT*, float*, KK_INT*);
void F77_BLAS_MANGLE(dsyr, DSYR)(const char*, KK_INT*, const double*, const double*, KK_INT*, double*, KK_INT*);
// Although there is a cgeru, there is no csyru
// Although there is a zgeru, there is no zsyru
// Although there is a cgerc, there is no csyrc, but there is cher (see below)
// Although there is a zgerc, there is no zsyrc, but there is zher (see below)

///
/// Her
///

void F77_BLAS_MANGLE(cher, CHER)(const char*, KK_INT*, const float*, const std::complex<float>*, KK_INT*,
                                 /* */ std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(zher, ZHER)(const char*, KK_INT*, const double*, const std::complex<double>*, KK_INT*,
                                 /* */ std::complex<double>*, KK_INT*);

///
/// Syr2
///
void F77_BLAS_MANGLE(ssyr2, SSYR2)(const char*, KK_INT*, const float*, const float*, const KK_INT*, const float*,
                                   KK_INT*, float*, KK_INT*);
void F77_BLAS_MANGLE(dsyr2, DSYR2)(const char*, KK_INT*, const double*, const double*, const KK_INT*, const double*,
                                   KK_INT*, double*, KK_INT*);
// Although there is a cgeru, there is no csyr2u
// Although there is a zgeru, there is no zsyr2u
// Although there is a cgerc, there is no csyr2c, but there is cher2 (see below)
// Although there is a zgerc, there is no zsyr2c, but there is zher2 (see below)

///
/// Her2
///

void F77_BLAS_MANGLE(cher2, CHER2)(const char*, KK_INT*, const std::complex<float>*, const std::complex<float>*,
                                   KK_INT*, const std::complex<float>*, KK_INT*,
                                   /* */ std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(zher2, ZHER2)(const char*, KK_INT*, const std::complex<double>*, const std::complex<double>*,
                                   KK_INT*, const std::complex<double>*, KK_INT*,
                                   /* */ std::complex<double>*, KK_INT*);

///
/// Trsv
///

void F77_BLAS_MANGLE(strsv, STRSV)(const char*, const char*, const char*, KK_INT*, const float*, KK_INT*,
                                   /* */ float*, KK_INT*);
void F77_BLAS_MANGLE(dtrsv, DTRSV)(const char*, const char*, const char*, KK_INT*, const double*, KK_INT*,
                                   /* */ double*, KK_INT*);
void F77_BLAS_MANGLE(ctrsv, CTRSV)(const char*, const char*, const char*, KK_INT*, const std::complex<float>*, KK_INT*,
                                   /* */ std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(ztrsv, ZTRSV)(const char*, const char*, const char*, KK_INT*, const std::complex<double>*, KK_INT*,
                                   /* */ std::complex<double>*, KK_INT*);

///
/// Gemm
///

void F77_BLAS_MANGLE(sgemm, SGEMM)(const char*, const char*, KK_INT*, KK_INT*, KK_INT*, const float*, const float*,
                                   KK_INT*, const float*, KK_INT*, const float*,
                                   /* */ float*, KK_INT*);
void F77_BLAS_MANGLE(dgemm, DGEMM)(const char*, const char*, KK_INT*, KK_INT*, KK_INT*, const double*, const double*,
                                   KK_INT*, const double*, KK_INT*, const double*,
                                   /* */ double*, KK_INT*);
void F77_BLAS_MANGLE(cgemm, CGEMM)(const char*, const char*, KK_INT*, KK_INT*, KK_INT*, const std::complex<float>*,
                                   const std::complex<float>*, KK_INT*, const std::complex<float>*, KK_INT*,
                                   const std::complex<float>*,
                                   /* */ std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(zgemm, ZGEMM)(const char*, const char*, KK_INT*, KK_INT*, KK_INT*, const std::complex<double>*,
                                   const std::complex<double>*, KK_INT*, const std::complex<double>*, KK_INT*,
                                   const std::complex<double>*,
                                   /* */ std::complex<double>*, KK_INT*);

///
/// Herk
///

void F77_BLAS_MANGLE(ssyrk, SSYRK)(const char*, const char*, KK_INT*, KK_INT*, const float*, const float*, KK_INT*,
                                   const float*,
                                   /* */ float*, KK_INT*);
void F77_BLAS_MANGLE(dsyrk, DSYRK)(const char*, const char*, KK_INT*, KK_INT*, const double*, const double*, KK_INT*,
                                   const double*,
                                   /* */ double*, KK_INT*);
void F77_BLAS_MANGLE(cherk, CHERK)(const char*, const char*, KK_INT*, KK_INT*, const std::complex<float>*,
                                   const std::complex<float>*, KK_INT*, const std::complex<float>*,
                                   /* */ std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(zherk, ZHERK)(const char*, const char*, KK_INT*, KK_INT*, const std::complex<double>*,
                                   const std::complex<double>*, KK_INT*, const std::complex<double>*,
                                   /* */ std::complex<double>*, KK_INT*);

///
/// Trmm
///

void F77_BLAS_MANGLE(strmm, STRMM)(const char*, const char*, const char*, const char*, KK_INT*, KK_INT*, const float*,
                                   const float*, KK_INT*,
                                   /* */ float*, KK_INT*);
void F77_BLAS_MANGLE(dtrmm, DTRMM)(const char*, const char*, const char*, const char*, KK_INT*, KK_INT*, const double*,
                                   const double*, KK_INT*,
                                   /* */ double*, KK_INT*);
void F77_BLAS_MANGLE(ctrmm, CTRMM)(const char*, const char*, const char*, const char*, KK_INT*, KK_INT*,
                                   const std::complex<float>*, const std::complex<float>*, KK_INT*,
                                   /* */ std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(ztrmm, ZTRMM)(const char*, const char*, const char*, const char*, KK_INT*, KK_INT*,
                                   const std::complex<double>*, const std::complex<double>*, KK_INT*,
                                   /* */ std::complex<double>*, KK_INT*);

///
/// Trsm
///

void F77_BLAS_MANGLE(strsm, STRSM)(const char*, const char*, const char*, const char*, KK_INT*, KK_INT*, const float*,
                                   const float*, KK_INT*,
                                   /* */ float*, KK_INT*);
void F77_BLAS_MANGLE(dtrsm, DTRSM)(const char*, const char*, const char*, const char*, KK_INT*, KK_INT*, const double*,
                                   const double*, KK_INT*,
                                   /* */ double*, KK_INT*);
void F77_BLAS_MANGLE(ctrsm, CTRSM)(const char*, const char*, const char*, const char*, KK_INT*, KK_INT*,
                                   const std::complex<float>*, const std::complex<float>*, KK_INT*,
                                   /* */ std::complex<float>*, KK_INT*);
void F77_BLAS_MANGLE(ztrsm, ZTRSM)(const char*, const char*, const char*, const char*, KK_INT*, KK_INT*,
                                   const std::complex<double>*, const std::complex<double>*, KK_INT*,
                                   /* */ std::complex<double>*, KK_INT*);

void F77_BLAS_MANGLE(sscal, SSCAL)(const KK_INT* N, const float* alpha,
                                   /* */ float* x, const KK_INT* x_inc);
void F77_BLAS_MANGLE(dscal, DSCAL)(const KK_INT* N, const double* alpha,
                                   /* */ double* x, const KK_INT* x_inc);
void F77_BLAS_MANGLE(cscal, CSCAL)(const KK_INT* N, const std::complex<float>* alpha,
                                   /* */ std::complex<float>* x, const KK_INT* x_inc);
void F77_BLAS_MANGLE(zscal, ZSCAL)(const KK_INT* N, const std::complex<double>* alpha,
                                   /* */ std::complex<double>* x, const KK_INT* x_inc);
}  // extern "C"

#define F77_FUNC_SSCAL F77_BLAS_MANGLE(sscal, SSCAL)
#define F77_FUNC_DSCAL F77_BLAS_MANGLE(dscal, DSCAL)
#define F77_FUNC_CSCAL F77_BLAS_MANGLE(cscal, CSCAL)
#define F77_FUNC_ZSCAL F77_BLAS_MANGLE(zscal, ZSCAL)

#define F77_FUNC_ISAMAX F77_BLAS_MANGLE(isamax, ISAMAX)
#define F77_FUNC_IDAMAX F77_BLAS_MANGLE(idamax, IDAMAX)
#define F77_FUNC_ICAMAX F77_BLAS_MANGLE(icamax, ICAMAX)
#define F77_FUNC_IZAMAX F77_BLAS_MANGLE(izamax, IZAMAX)

#define F77_FUNC_SNRM2 F77_BLAS_MANGLE(snrm2, SNRM2)
#define F77_FUNC_DNRM2 F77_BLAS_MANGLE(dnrm2, DNRM2)
#define F77_FUNC_SCNRM2 F77_BLAS_MANGLE(scnrm2, SCNRM2)
#define F77_FUNC_DZNRM2 F77_BLAS_MANGLE(dznrm2, DZNRM2)

#define F77_FUNC_SASUM F77_BLAS_MANGLE(sasum, SASUM)
#define F77_FUNC_DASUM F77_BLAS_MANGLE(dasum, DASUM)
#define F77_FUNC_SCASUM F77_BLAS_MANGLE(scasum, SCASUM)
#define F77_FUNC_DZASUM F77_BLAS_MANGLE(dzasum, DZASUM)

#define F77_FUNC_SDOT F77_BLAS_MANGLE(sdot, SDOT)
#define F77_FUNC_DDOT F77_BLAS_MANGLE(ddot, DDOT)
#define F77_FUNC_CDOTU F77_BLAS_MANGLE(cdotu, CDOTU)
#define F77_FUNC_ZDOTU F77_BLAS_MANGLE(zdotu, ZDOTU)
#define F77_FUNC_CDOTC F77_BLAS_MANGLE(cdotc, CDOTC)
#define F77_FUNC_ZDOTC F77_BLAS_MANGLE(zdotc, ZDOTC)

#define F77_FUNC_SAXPY F77_BLAS_MANGLE(saxpy, SAXPY)
#define F77_FUNC_DAXPY F77_BLAS_MANGLE(daxpy, DAXPY)
#define F77_FUNC_CAXPY F77_BLAS_MANGLE(caxpy, CAXPY)
#define F77_FUNC_ZAXPY F77_BLAS_MANGLE(zaxpy, ZAXPY)

#define F77_FUNC_SROT F77_BLAS_MANGLE(srot, SROT)
#define F77_FUNC_DROT F77_BLAS_MANGLE(drot, DROT)
#define F77_FUNC_CROT F77_BLAS_MANGLE(crot, CROT)
#define F77_FUNC_ZROT F77_BLAS_MANGLE(zrot, ZROT)

#define F77_FUNC_SROTG F77_BLAS_MANGLE(srotg, SROTG)
#define F77_FUNC_DROTG F77_BLAS_MANGLE(drotg, DROTG)
#define F77_FUNC_CROTG F77_BLAS_MANGLE(crotg, CROTG)
#define F77_FUNC_ZROTG F77_BLAS_MANGLE(zrotg, ZROTG)

#define F77_FUNC_SROTM F77_BLAS_MANGLE(srotm, SROTM)
#define F77_FUNC_DROTM F77_BLAS_MANGLE(drotm, DROTM)

#define F77_FUNC_SROTMG F77_BLAS_MANGLE(srotmg, SROTMG)
#define F77_FUNC_DROTMG F77_BLAS_MANGLE(drotmg, DROTMG)

#define F77_FUNC_SSWAP F77_BLAS_MANGLE(sswap, SSWAP)
#define F77_FUNC_DSWAP F77_BLAS_MANGLE(dswap, DSWAP)
#define F77_FUNC_CSWAP F77_BLAS_MANGLE(cswap, CSWAP)
#define F77_FUNC_ZSWAP F77_BLAS_MANGLE(zswap, ZSWAP)

#define F77_FUNC_SGEMV F77_BLAS_MANGLE(sgemv, SGEMV)
#define F77_FUNC_DGEMV F77_BLAS_MANGLE(dgemv, DGEMV)
#define F77_FUNC_CGEMV F77_BLAS_MANGLE(cgemv, CGEMV)
#define F77_FUNC_ZGEMV F77_BLAS_MANGLE(zgemv, ZGEMV)

#define F77_FUNC_SGER F77_BLAS_MANGLE(sger, SGER)
#define F77_FUNC_DGER F77_BLAS_MANGLE(dger, DGER)
#define F77_FUNC_CGERU F77_BLAS_MANGLE(cgeru, CGERU)
#define F77_FUNC_ZGERU F77_BLAS_MANGLE(zgeru, ZGERU)
#define F77_FUNC_CGERC F77_BLAS_MANGLE(cgerc, CGERC)
#define F77_FUNC_ZGERC F77_BLAS_MANGLE(zgerc, ZGERC)

#define F77_FUNC_SSYR F77_BLAS_MANGLE(ssyr, SSYR)
#define F77_FUNC_DSYR F77_BLAS_MANGLE(dsyr, DSYR)

#define F77_FUNC_CHER F77_BLAS_MANGLE(cher, CHER)
#define F77_FUNC_ZHER F77_BLAS_MANGLE(zher, ZHER)

#define F77_FUNC_SSYR2 F77_BLAS_MANGLE(ssyr2, SSYR2)
#define F77_FUNC_DSYR2 F77_BLAS_MANGLE(dsyr2, DSYR2)

#define F77_FUNC_CHER2 F77_BLAS_MANGLE(cher2, CHER2)
#define F77_FUNC_ZHER2 F77_BLAS_MANGLE(zher2, ZHER2)

#define F77_FUNC_STRSV F77_BLAS_MANGLE(strsv, STRSV)
#define F77_FUNC_DTRSV F77_BLAS_MANGLE(dtrsv, DTRSV)
#define F77_FUNC_CTRSV F77_BLAS_MANGLE(ctrsv, CTRSV)
#define F77_FUNC_ZTRSV F77_BLAS_MANGLE(ztrsv, ZTRSV)

#define F77_FUNC_SGEMM F77_BLAS_MANGLE(sgemm, SGEMM)
#define F77_FUNC_DGEMM F77_BLAS_MANGLE(dgemm, DGEMM)
#define F77_FUNC_CGEMM F77_BLAS_MANGLE(cgemm, CGEMM)
#define F77_FUNC_ZGEMM F77_BLAS_MANGLE(zgemm, ZGEMM)

#define F77_FUNC_SSYRK F77_BLAS_MANGLE(ssyrk, SSYRK)
#define F77_FUNC_DSYRK F77_BLAS_MANGLE(dsyrk, DSYRK)
#define F77_FUNC_CHERK F77_BLAS_MANGLE(cherk, CHERK)
#define F77_FUNC_ZHERK F77_BLAS_MANGLE(zherk, ZHERK)

#define F77_FUNC_STRMM F77_BLAS_MANGLE(strmm, STRMM)
#define F77_FUNC_DTRMM F77_BLAS_MANGLE(dtrmm, DTRMM)
#define F77_FUNC_CTRMM F77_BLAS_MANGLE(ctrmm, CTRMM)
#define F77_FUNC_ZTRMM F77_BLAS_MANGLE(ztrmm, ZTRMM)

#define F77_FUNC_STRSM F77_BLAS_MANGLE(strsm, STRSM)
#define F77_FUNC_DTRSM F77_BLAS_MANGLE(dtrsm, DTRSM)
#define F77_FUNC_CTRSM F77_BLAS_MANGLE(ctrsm, CTRSM)
#define F77_FUNC_ZTRSM F77_BLAS_MANGLE(ztrsm, ZTRSM)

#endif  // defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)

namespace KokkosBlas {
namespace Impl {

///
/// float
///

template <>
void HostBlas<float>::scal(KK_INT n, const float alpha,
                           /* */ float* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_sscal(n, alpha, x, x_inc);
#else
  F77_FUNC_SSCAL(&n, &alpha, x, &x_inc);
#endif
}
template <>
KK_INT HostBlas<float>::iamax(KK_INT n, const float* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_isamax(n, x, x_inc) + 1;
#else
  return F77_FUNC_ISAMAX(&n, x, &x_inc);
#endif
}
template <>
float HostBlas<float>::nrm2(KK_INT n, const float* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_snrm2(n, x, x_inc);
#else
  return F77_FUNC_SNRM2(&n, x, &x_inc);
#endif
}
template <>
float HostBlas<float>::asum(KK_INT n, const float* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_sasum(n, x, x_inc);
#else
  return F77_FUNC_SASUM(&n, x, &x_inc);
#endif
}
template <>
float HostBlas<float>::dot(KK_INT n, const float* x, KK_INT x_inc, const float* y, KK_INT y_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_sdot(n, x, x_inc, y, y_inc);
#else
  return F77_FUNC_SDOT(&n, x, &x_inc, y, &y_inc);
#endif
}
template <>
void HostBlas<float>::axpy(KK_INT n, const float alpha, const float* x, KK_INT x_inc,
                           /* */ float* y, KK_INT y_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_saxpy(n, alpha, x, x_inc, y, y_inc);
#else
  F77_FUNC_SAXPY(&n, &alpha, x, &x_inc, y, &y_inc);
#endif
}
template <>
void HostBlas<float>::rot(KK_INT const N, float* X, KK_INT const incx, float* Y, KK_INT const incy, float* c,
                          float* s) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_srot(N, X, incx, Y, incy, *c, *s);
#else
  F77_FUNC_SROT(&N, X, &incx, Y, &incy, c, s);
#endif
}
template <>
void HostBlas<float>::rotg(float* a, float* b, float* c, float* s) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_srotg(a, b, c, s);
#else
  F77_FUNC_SROTG(a, b, c, s);
#endif
}
template <>
void HostBlas<float>::rotm(const KK_INT n, float* X, const KK_INT incx, float* Y, const KK_INT incy,
                           const float* param) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_srotm(n, X, incx, Y, incy, param);
#else
  F77_FUNC_SROTM(&n, X, &incx, Y, &incy, param);
#endif
}
template <>
void HostBlas<float>::rotmg(float* d1, float* d2, float* x1, const float* y1, float* param) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_srotmg(d1, d2, x1, *y1, param);
#else
  F77_FUNC_SROTMG(d1, d2, x1, y1, param);
#endif
}
template <>
void HostBlas<float>::swap(KK_INT const N, float* X, KK_INT const incx, float* Y, KK_INT const incy) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_sswap(N, X, incx, Y, incy);
#else
  F77_FUNC_SSWAP(&N, X, &incx, Y, &incy);
#endif
}
template <>
void HostBlas<float>::gemv(const char trans, KK_INT m, KK_INT n, const float alpha, const float* a, KK_INT lda,
                           const float* x, KK_INT incx, const float beta,
                           /* */ float* y, KK_INT incy) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order     = CblasColMajor;
  enum CBLAS_TRANSPOSE acc_trans = kk_to_accelerate_trans(trans);
  cblas_sgemv(acc_order, acc_trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
  F77_FUNC_SGEMV(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endif
}
template <>
void HostBlas<float>::ger(KK_INT m, KK_INT n, const float alpha, const float* x, KK_INT incx, const float* y,
                          KK_INT incy, float* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  cblas_sger(acc_order, m, n, alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_SGER(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
#endif
}
template <>
void HostBlas<float>::syr(const char uplo, KK_INT n, const float alpha, const float* x, KK_INT incx, float* a,
                          KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  cblas_ssyr(acc_order, acc_uplo, n, alpha, x, incx, a, lda);
#else
  F77_FUNC_SSYR(&uplo, &n, &alpha, x, &incx, a, &lda);
#endif
}
template <>
void HostBlas<float>::syr2(const char uplo, KK_INT n, const float alpha, const float* x, KK_INT incx, const float* y,
                           KK_INT incy, float* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  cblas_ssyr2(acc_order, acc_uplo, n, alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_SSYR2(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
#endif
}
template <>
void HostBlas<float>::trsv(const char uplo, const char transa, const char diag, KK_INT m, const float* a, KK_INT lda,
                           /* */ float* x, KK_INT incx) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_strsv(acc_order, acc_uplo, acc_trans, acc_diag, m, a, lda, x, incx);
#else
  F77_FUNC_STRSV(&uplo, &transa, &diag, &m, a, &lda, x, &incx);
#endif
}
template <>
void HostBlas<float>::gemm(const char transa, const char transb, KK_INT m, KK_INT n, KK_INT k, const float alpha,
                           const float* a, KK_INT lda, const float* b, KK_INT ldb, const float beta,
                           /* */ float* c, KK_INT ldc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  auto acc_transa            = kk_to_accelerate_trans(transa);
  auto acc_transb            = kk_to_accelerate_trans(transb);
  cblas_sgemm(acc_order, acc_transa, acc_transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  F77_FUNC_SGEMM(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
}
template <>
void HostBlas<float>::herk(const char uplo, const char trans, KK_INT n, KK_INT k, const float alpha, const float* a,
                           KK_INT lda, const float beta,
                           /* */ float* c, KK_INT ldc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(trans);
  cblas_ssyrk(acc_order, acc_uplo, acc_trans, n, k, alpha, a, lda, beta, c, ldc);
#else
  F77_FUNC_SSYRK(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
#endif
}
template <>
void HostBlas<float>::trmm(const char side, const char uplo, const char transa, const char diag, KK_INT m, KK_INT n,
                           const float alpha, const float* a, KK_INT lda,
                           /* */ float* b, KK_INT ldb) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  CBLAS_SIDE acc_side        = kk_to_accelerate_side(side);
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_strmm(acc_order, acc_side, acc_uplo, acc_trans, acc_diag, m, n, alpha, a, lda, b, ldb);
#else
  F77_FUNC_STRMM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
#endif
}
template <>
void HostBlas<float>::trsm(const char side, const char uplo, const char transa, const char diag, KK_INT m, KK_INT n,
                           const float alpha, const float* a, KK_INT lda,
                           /* */ float* b, KK_INT ldb) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  CBLAS_SIDE acc_side        = kk_to_accelerate_side(side);
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_strsm(acc_order, acc_side, acc_uplo, acc_trans, acc_diag, m, n, alpha, a, lda, b, ldb);
#else
  F77_FUNC_STRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
#endif
}

///
/// double
///

template <>
void HostBlas<double>::scal(KK_INT n, const double alpha,
                            /* */ double* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_dscal(n, alpha, x, x_inc);
#else
  F77_FUNC_DSCAL(&n, &alpha, x, &x_inc);
#endif
}
template <>
KK_INT HostBlas<double>::iamax(KK_INT n, const double* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_idamax(n, x, x_inc) + 1;
#else
  return F77_FUNC_IDAMAX(&n, x, &x_inc);
#endif
}
template <>
double HostBlas<double>::nrm2(KK_INT n, const double* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_dnrm2(n, x, x_inc);
#else
  return F77_FUNC_DNRM2(&n, x, &x_inc);
#endif
}
template <>
double HostBlas<double>::asum(KK_INT n, const double* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_dasum(n, x, x_inc);
#else
  return F77_FUNC_DASUM(&n, x, &x_inc);
#endif
}
template <>
double HostBlas<double>::dot(KK_INT n, const double* x, KK_INT x_inc, const double* y, KK_INT y_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_ddot(n, x, x_inc, y, y_inc);
#else
  return F77_FUNC_DDOT(&n, x, &x_inc, y, &y_inc);
#endif
}
template <>
void HostBlas<double>::axpy(KK_INT n, const double alpha, const double* x, KK_INT x_inc,
                            /* */ double* y, KK_INT y_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_daxpy(n, alpha, x, x_inc, y, y_inc);
#else
  F77_FUNC_DAXPY(&n, &alpha, x, &x_inc, y, &y_inc);
#endif
}
template <>
void HostBlas<double>::rot(KK_INT const N, double* X, KK_INT const incx, double* Y, KK_INT const incy, double* c,
                           double* s) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_drot(N, X, incx, Y, incy, *c, *s);
#else
  F77_FUNC_DROT(&N, X, &incx, Y, &incy, c, s);
#endif
}
template <>
void HostBlas<double>::rotg(double* a, double* b, double* c, double* s) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_drotg(a, b, c, s);
#else
  F77_FUNC_DROTG(a, b, c, s);
#endif
}
template <>
void HostBlas<double>::rotm(const KK_INT n, double* X, const KK_INT incx, double* Y, const KK_INT incy,
                            const double* param) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_drotm(n, X, incx, Y, incy, param);
#else
  F77_FUNC_DROTM(&n, X, &incx, Y, &incy, param);
#endif
}
template <>
void HostBlas<double>::rotmg(double* d1, double* d2, double* x1, const double* y1, double* param) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_drotmg(d1, d2, x1, *y1, param);
#else
  F77_FUNC_DROTMG(d1, d2, x1, y1, param);
#endif
}
template <>
void HostBlas<double>::swap(KK_INT const N, double* X, KK_INT const incx, double* Y, KK_INT const incy) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_dswap(N, X, incx, Y, incy);
#else
  F77_FUNC_DSWAP(&N, X, &incx, Y, &incy);
#endif
}
template <>
void HostBlas<double>::gemv(const char trans, KK_INT m, KK_INT n, const double alpha, const double* a, KK_INT lda,
                            const double* x, KK_INT incx, const double beta,
                            /* */ double* y, KK_INT incy) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order     = CblasColMajor;
  enum CBLAS_TRANSPOSE acc_trans = kk_to_accelerate_trans(trans);
  cblas_dgemv(acc_order, acc_trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
  F77_FUNC_DGEMV(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endif
}
template <>
void HostBlas<double>::ger(KK_INT m, KK_INT n, const double alpha, const double* x, KK_INT incx, const double* y,
                           KK_INT incy, double* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  cblas_dger(acc_order, m, n, alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_DGER(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
#endif
}
template <>
void HostBlas<double>::syr(const char uplo, KK_INT n, const double alpha, const double* x, KK_INT incx, double* a,
                           KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  cblas_dsyr(acc_order, acc_uplo, n, alpha, x, incx, a, lda);
#else
  F77_FUNC_DSYR(&uplo, &n, &alpha, x, &incx, a, &lda);
#endif
}
template <>
void HostBlas<double>::syr2(const char uplo, KK_INT n, const double alpha, const double* x, KK_INT incx,
                            const double* y, KK_INT incy, double* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  cblas_dsyr2(acc_order, acc_uplo, n, alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_DSYR2(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
#endif
}
template <>
void HostBlas<double>::trsv(const char uplo, const char transa, const char diag, KK_INT m, const double* a, KK_INT lda,
                            /* */ double* x, KK_INT incx) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_dtrsv(acc_order, acc_uplo, acc_trans, acc_diag, m, a, lda, x, incx);
#else
  F77_FUNC_DTRSV(&uplo, &transa, &diag, &m, a, &lda, x, &incx);
#endif
}
template <>
void HostBlas<double>::gemm(const char transa, const char transb, KK_INT m, KK_INT n, KK_INT k, const double alpha,
                            const double* a, KK_INT lda, const double* b, KK_INT ldb, const double beta,
                            /* */ double* c, KK_INT ldc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  auto acc_transa            = kk_to_accelerate_trans(transa);
  auto acc_transb            = kk_to_accelerate_trans(transb);
  cblas_dgemm(acc_order, acc_transa, acc_transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  F77_FUNC_DGEMM(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
}
template <>
void HostBlas<double>::herk(const char uplo, const char trans, KK_INT n, KK_INT k, const double alpha, const double* a,
                            KK_INT lda, const double beta,
                            /* */ double* c, KK_INT ldc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(trans);
  cblas_dsyrk(acc_order, acc_uplo, acc_trans, n, k, alpha, a, lda, beta, c, ldc);
#else
  F77_FUNC_DSYRK(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
#endif
}
template <>
void HostBlas<double>::trmm(const char side, const char uplo, const char transa, const char diag, KK_INT m, KK_INT n,
                            const double alpha, const double* a, KK_INT lda,
                            /* */ double* b, KK_INT ldb) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  CBLAS_SIDE acc_side        = kk_to_accelerate_side(side);
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_dtrmm(acc_order, acc_side, acc_uplo, acc_trans, acc_diag, m, n, alpha, a, lda, b, ldb);
#else
  F77_FUNC_DTRMM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
#endif
}
template <>
void HostBlas<double>::trsm(const char side, const char uplo, const char transa, const char diag, KK_INT m, KK_INT n,
                            const double alpha, const double* a, KK_INT lda,
                            /* */ double* b, KK_INT ldb) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  CBLAS_SIDE acc_side        = kk_to_accelerate_side(side);
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_dtrsm(acc_order, acc_side, acc_uplo, acc_trans, acc_diag, m, n, alpha, a, lda, b, ldb);
#else
  F77_FUNC_DTRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
#endif
}

///
/// std::complex<float>
///

template <>
void HostBlas<std::complex<float> >::scal(KK_INT n, const std::complex<float> alpha,
                                          /* */ std::complex<float>* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_cscal(n, &alpha, x, x_inc);
#else
  F77_FUNC_CSCAL(&n, &alpha, x, &x_inc);
#endif
}
template <>
KK_INT HostBlas<std::complex<float> >::iamax(KK_INT n, const std::complex<float>* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_icamax(n, x, x_inc) + 1;
#else
  return F77_FUNC_ICAMAX(&n, x, &x_inc);
#endif
}
template <>
float HostBlas<std::complex<float> >::nrm2(KK_INT n, const std::complex<float>* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_scnrm2(n, x, x_inc);
#else
  return F77_FUNC_SCNRM2(&n, x, &x_inc);
#endif
}
template <>
float HostBlas<std::complex<float> >::asum(KK_INT n, const std::complex<float>* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_scasum(n, x, x_inc);
#else
  return F77_FUNC_SCASUM(&n, x, &x_inc);
#endif
}
template <>
std::complex<float> HostBlas<std::complex<float> >::dot(KK_INT n, const std::complex<float>* x, KK_INT x_inc,
                                                        const std::complex<float>* y, KK_INT y_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  std::complex<float> res;
  cblas_cdotc_sub(n, x, x_inc, y, y_inc, &res);
  return res;
#else
#if defined(KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX)
  _kk_float2 res = F77_FUNC_CDOTC(&n, x, &x_inc, y, &y_inc);
  return std::complex<float>(res.vals[0], res.vals[1]);
#else
  std::complex<float> res;
  F77_FUNC_CDOTC(&res, &n, x, &x_inc, y, &y_inc);
  return res;
#endif
#endif
}
template <>
void HostBlas<std::complex<float> >::axpy(KK_INT n, const std::complex<float> alpha, const std::complex<float>* x,
                                          KK_INT x_inc,
                                          /* */ std::complex<float>* y, KK_INT y_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_caxpy(n, &alpha, x, x_inc, y, y_inc);
#else
  F77_FUNC_CAXPY(&n, &alpha, x, &x_inc, y, &y_inc);
#endif
}
template <>
void HostBlas<std::complex<float> >::rot(KK_INT const N, std::complex<float>* X, KK_INT const incx,
                                         std::complex<float>* Y, KK_INT const incy, float* c, std::complex<float>* s) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_csrot(N, X, incx, Y, incy, *c, s->real());
#else
  F77_FUNC_CROT(&N, X, &incx, Y, &incy, c, s);
#endif
}
template <>
void HostBlas<std::complex<float> >::rotg(std::complex<float>* a, std::complex<float>* b, float* c,
                                          std::complex<float>* s) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_crotg(a, b, c, s);
#else
  F77_FUNC_CROTG(a, b, c, s);
#endif
}
template <>
void HostBlas<std::complex<float> >::swap(KK_INT const N, std::complex<float>* X, KK_INT const incx,
                                          std::complex<float>* Y, KK_INT const incy) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_cswap(N, X, incx, Y, incy);
#else
  F77_FUNC_CSWAP(&N, X, &incx, Y, &incy);
#endif
}
template <>
void HostBlas<std::complex<float> >::gemv(const char trans, KK_INT m, KK_INT n, const std::complex<float> alpha,
                                          const std::complex<float>* a, KK_INT lda, const std::complex<float>* x,
                                          KK_INT incx, const std::complex<float> beta,
                                          /* */ std::complex<float>* y, KK_INT incy) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order     = CblasColMajor;
  enum CBLAS_TRANSPOSE acc_trans = kk_to_accelerate_trans(trans);
  cblas_cgemv(acc_order, acc_trans, m, n, &alpha, a, lda, x, incx, &beta, y, incy);
#else
  F77_FUNC_CGEMV(&trans, &m, &n, &alpha, (const std::complex<float>*)a, &lda, (const std::complex<float>*)x, &incx,
                 &beta, (std::complex<float>*)y, &incy);
#endif
}
template <>
void HostBlas<std::complex<float> >::geru(KK_INT m, KK_INT n, const std::complex<float> alpha,
                                          const std::complex<float>* x, KK_INT incx, const std::complex<float>* y,
                                          KK_INT incy, std::complex<float>* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  cblas_cgeru(acc_order, m, n, &alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_CGERU(&m, &n, &alpha, (const std::complex<float>*)x, &incx, (const std::complex<float>*)y, &incy,
                 (std::complex<float>*)a, &lda);
#endif
}
template <>
void HostBlas<std::complex<float> >::gerc(KK_INT m, KK_INT n, const std::complex<float> alpha,
                                          const std::complex<float>* x, KK_INT incx, const std::complex<float>* y,
                                          KK_INT incy, std::complex<float>* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  cblas_cgerc(acc_order, m, n, &alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_CGERC(&m, &n, &alpha, (const std::complex<float>*)x, &incx, (const std::complex<float>*)y, &incy,
                 (std::complex<float>*)a, &lda);
#endif
}
template <>
template <>
void HostBlas<std::complex<float> >::her<float>(const char uplo, KK_INT n, const float alpha,
                                                const std::complex<float>* x, KK_INT incx, std::complex<float>* a,
                                                KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  cblas_cher(acc_order, acc_uplo, n, alpha, x, incx, a, lda);
#else
  F77_FUNC_CHER(&uplo, &n, &alpha, (const std::complex<float>*)x, &incx, (std::complex<float>*)a, &lda);
#endif
}
template <>
void HostBlas<std::complex<float> >::her2(const char uplo, KK_INT n, const std::complex<float> alpha,
                                          const std::complex<float>* x, KK_INT incx, const std::complex<float>* y,
                                          KK_INT incy, std::complex<float>* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  cblas_cher2(acc_order, acc_uplo, n, &alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_CHER2(&uplo, &n, &alpha, (const std::complex<float>*)x, &incx, (const std::complex<float>*)y, &incy,
                 (std::complex<float>*)a, &lda);
#endif
}
template <>
void HostBlas<std::complex<float> >::trsv(const char uplo, const char transa, const char diag, KK_INT m,
                                          const std::complex<float>* a, KK_INT lda,
                                          /* */ std::complex<float>* x, KK_INT incx) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_ctrsv(acc_order, acc_uplo, acc_trans, acc_diag, m, a, lda, x, incx);
#else
  F77_FUNC_CTRSV(&uplo, &transa, &diag, &m, (const std::complex<float>*)a, &lda, (std::complex<float>*)x, &incx);
#endif
}
template <>
void HostBlas<std::complex<float> >::gemm(const char transa, const char transb, KK_INT m, KK_INT n, KK_INT k,
                                          const std::complex<float> alpha, const std::complex<float>* a, KK_INT lda,
                                          const std::complex<float>* b, KK_INT ldb, const std::complex<float> beta,
                                          /* */ std::complex<float>* c, KK_INT ldc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  auto acc_transa            = kk_to_accelerate_trans(transa);
  auto acc_transb            = kk_to_accelerate_trans(transb);
  cblas_cgemm(acc_order, acc_transa, acc_transb, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
#else
  F77_FUNC_CGEMM(&transa, &transb, &m, &n, &k, &alpha, (const std::complex<float>*)a, &lda,
                 (const std::complex<float>*)b, &ldb, &beta, (std::complex<float>*)c, &ldc);
#endif
}
template <>
void HostBlas<std::complex<float> >::herk(const char uplo, const char transa, KK_INT n, KK_INT k,
                                          const std::complex<float> alpha, const std::complex<float>* a, KK_INT lda,
                                          const std::complex<float> beta,
                                          /* */ std::complex<float>* c, KK_INT ldc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  cblas_cherk(acc_order, acc_uplo, acc_trans, n, k, alpha.real(), a, lda, beta.real(), c, ldc);
#else
  F77_FUNC_CHERK(&uplo, &transa, &n, &k, &alpha, (const std::complex<float>*)a, &lda, &beta, (std::complex<float>*)c,
                 &ldc);
#endif
}
template <>
void HostBlas<std::complex<float> >::trmm(const char side, const char uplo, const char transa, const char diag,
                                          KK_INT m, KK_INT n, const std::complex<float> alpha,
                                          const std::complex<float>* a, KK_INT lda,
                                          /* */ std::complex<float>* b, KK_INT ldb) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  CBLAS_SIDE acc_side        = kk_to_accelerate_side(side);
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_ctrmm(acc_order, acc_side, acc_uplo, acc_trans, acc_diag, m, n, &alpha, a, lda, b, ldb);
#else
  F77_FUNC_CTRMM(&side, &uplo, &transa, &diag, &m, &n, &alpha, (const std::complex<float>*)a, &lda,
                 (std::complex<float>*)b, &ldb);
#endif
}
template <>
void HostBlas<std::complex<float> >::trsm(const char side, const char uplo, const char transa, const char diag,
                                          KK_INT m, KK_INT n, const std::complex<float> alpha,
                                          const std::complex<float>* a, KK_INT lda,
                                          /* */ std::complex<float>* b, KK_INT ldb) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  CBLAS_SIDE acc_side        = kk_to_accelerate_side(side);
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_ctrsm(acc_order, acc_side, acc_uplo, acc_trans, acc_diag, m, n, &alpha, a, lda, b, ldb);
#else
  F77_FUNC_CTRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, (const std::complex<float>*)a, &lda,
                 (std::complex<float>*)b, &ldb);
#endif
}

///
/// std::complex<double>
///

template <>
void HostBlas<std::complex<double> >::scal(KK_INT n, const std::complex<double> alpha,
                                           /* */ std::complex<double>* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_zscal(n, &alpha, x, x_inc);
#else
  F77_FUNC_ZSCAL(&n, &alpha, x, &x_inc);
#endif
}
template <>
KK_INT HostBlas<std::complex<double> >::iamax(KK_INT n, const std::complex<double>* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_izamax(n, x, x_inc) + 1;
#else
  return F77_FUNC_IZAMAX(&n, x, &x_inc);
#endif
}
template <>
double HostBlas<std::complex<double> >::nrm2(KK_INT n, const std::complex<double>* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_dznrm2(n, x, x_inc);
#else
  return F77_FUNC_DZNRM2(&n, x, &x_inc);
#endif
}
template <>
double HostBlas<std::complex<double> >::asum(KK_INT n, const std::complex<double>* x, KK_INT x_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  return cblas_dzasum(n, x, x_inc);
#else
  // see issue 2005
  // On some platforms with OpenBLAS < 0.3.26, dzasum on vectors less than 16 entries is producing 0.
  // this has been observed on some (not all) systems with:
  // clang 14.0.6 / 15.0.7 AND OpenBLAS 0.3.23 AND Sapphire Rapids CPU
  // unfortunately, it's not clear exactly what the trigger is
  if (n > 0 && n < 16) {
    double ret = 0.0;
    for (int i = 0; i < n; ++i) {
      ret += Kokkos::abs(x[i].real()) + Kokkos::abs(x[i].imag());
    }
    return ret;
  }
  return F77_FUNC_DZASUM(&n, x, &x_inc);
#endif
}
template <>
std::complex<double> HostBlas<std::complex<double> >::dot(KK_INT n, const std::complex<double>* x, KK_INT x_inc,
                                                          const std::complex<double>* y, KK_INT y_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  std::complex<double> res;
  cblas_zdotc_sub(n, x, x_inc, y, y_inc, &res);
  return res;
#else
#if defined(KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX)
  _kk_double2 res = F77_FUNC_ZDOTC(&n, x, &x_inc, y, &y_inc);
  return std::complex<double>(res.vals[0], res.vals[1]);
#else
  std::complex<double> res;
  F77_FUNC_ZDOTC(&res, &n, x, &x_inc, y, &y_inc);
  return res;
#endif
#endif
}
template <>
void HostBlas<std::complex<double> >::axpy(KK_INT n, const std::complex<double> alpha, const std::complex<double>* x,
                                           KK_INT x_inc,
                                           /* */ std::complex<double>* y, KK_INT y_inc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_zaxpy(n, &alpha, x, x_inc, y, y_inc);
#else
  F77_FUNC_ZAXPY(&n, &alpha, x, &x_inc, y, &y_inc);
#endif
}
template <>
void HostBlas<std::complex<double> >::rot(KK_INT const N, std::complex<double>* X, KK_INT const incx,
                                          std::complex<double>* Y, KK_INT const incy, double* c,
                                          std::complex<double>* s) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_zdrot(N, X, incx, Y, incy, *c, s->real());
#else
  F77_FUNC_ZROT(&N, X, &incx, Y, &incy, c, s);
#endif
}
template <>
void HostBlas<std::complex<double> >::rotg(std::complex<double>* a, std::complex<double>* b, double* c,
                                           std::complex<double>* s) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_zrotg(a, b, c, s);
#else
  F77_FUNC_ZROTG(a, b, c, s);
#endif
}
template <>
void HostBlas<std::complex<double> >::swap(KK_INT const N, std::complex<double>* X, KK_INT const incx,
                                           std::complex<double>* Y, KK_INT const incy) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  cblas_zswap(N, X, incx, Y, incy);
#else
  F77_FUNC_ZSWAP(&N, X, &incx, Y, &incy);
#endif
}
template <>
void HostBlas<std::complex<double> >::gemv(const char trans, KK_INT m, KK_INT n, const std::complex<double> alpha,
                                           const std::complex<double>* a, KK_INT lda, const std::complex<double>* x,
                                           KK_INT incx, const std::complex<double> beta,
                                           /* */ std::complex<double>* y, KK_INT incy) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order     = CblasColMajor;
  enum CBLAS_TRANSPOSE acc_trans = kk_to_accelerate_trans(trans);
  cblas_zgemv(acc_order, acc_trans, m, n, &alpha, a, lda, x, incx, &beta, y, incy);
#else
  F77_FUNC_ZGEMV(&trans, &m, &n, &alpha, (const std::complex<double>*)a, &lda, (const std::complex<double>*)x, &incx,
                 &beta, (std::complex<double>*)y, &incy);
#endif
}
template <>
void HostBlas<std::complex<double> >::geru(KK_INT m, KK_INT n, const std::complex<double> alpha,
                                           const std::complex<double>* x, KK_INT incx, const std::complex<double>* y,
                                           KK_INT incy, std::complex<double>* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  cblas_zgeru(acc_order, m, n, &alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_ZGERU(&m, &n, &alpha, (const std::complex<double>*)x, &incx, (const std::complex<double>*)y, &incy,
                 (std::complex<double>*)a, &lda);
#endif
}
template <>
void HostBlas<std::complex<double> >::gerc(KK_INT m, KK_INT n, const std::complex<double> alpha,
                                           const std::complex<double>* x, KK_INT incx, const std::complex<double>* y,
                                           KK_INT incy, std::complex<double>* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  cblas_zgerc(acc_order, m, n, &alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_ZGERC(&m, &n, &alpha, (const std::complex<double>*)x, &incx, (const std::complex<double>*)y, &incy,
                 (std::complex<double>*)a, &lda);
#endif
}
template <>
template <>
void HostBlas<std::complex<double> >::her<double>(const char uplo, KK_INT n, const double alpha,
                                                  const std::complex<double>* x, KK_INT incx, std::complex<double>* a,
                                                  KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  cblas_zher(acc_order, acc_uplo, n, alpha, x, incx, a, lda);
#else
  F77_FUNC_ZHER(&uplo, &n, &alpha, (const std::complex<double>*)x, &incx, (std::complex<double>*)a, &lda);
#endif
}
template <>
void HostBlas<std::complex<double> >::her2(const char uplo, KK_INT n, const std::complex<double> alpha,
                                           const std::complex<double>* x, KK_INT incx, const std::complex<double>* y,
                                           KK_INT incy, std::complex<double>* a, KK_INT lda) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  cblas_zher2(acc_order, acc_uplo, n, &alpha, x, incx, y, incy, a, lda);
#else
  F77_FUNC_ZHER2(&uplo, &n, &alpha, (const std::complex<double>*)x, &incx, (const std::complex<double>*)y, &incy,
                 (std::complex<double>*)a, &lda);
#endif
}
template <>
void HostBlas<std::complex<double> >::trsv(const char uplo, const char transa, const char diag, KK_INT m,
                                           const std::complex<double>* a, KK_INT lda,
                                           /* */ std::complex<double>* x, KK_INT incx) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_ztrsv(acc_order, acc_uplo, acc_trans, acc_diag, m, a, lda, x, incx);
#else
  F77_FUNC_ZTRSV(&uplo, &transa, &diag, &m, (const std::complex<double>*)a, &lda, (std::complex<double>*)x, &incx);
#endif
}

template <>
void HostBlas<std::complex<double> >::gemm(const char transa, const char transb, KK_INT m, KK_INT n, KK_INT k,
                                           const std::complex<double> alpha, const std::complex<double>* a, KK_INT lda,
                                           const std::complex<double>* b, KK_INT ldb, const std::complex<double> beta,
                                           /* */ std::complex<double>* c, KK_INT ldc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  auto acc_transa            = kk_to_accelerate_trans(transa);
  auto acc_transb            = kk_to_accelerate_trans(transb);
  cblas_zgemm(acc_order, acc_transa, acc_transb, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
#else
  F77_FUNC_ZGEMM(&transa, &transb, &m, &n, &k, &alpha, (const std::complex<double>*)a, &lda,
                 (const std::complex<double>*)b, &ldb, &beta, (std::complex<double>*)c, &ldc);
#endif
}
template <>
void HostBlas<std::complex<double> >::herk(const char uplo, const char transa, KK_INT n, KK_INT k,
                                           const std::complex<double> alpha, const std::complex<double>* a, KK_INT lda,
                                           const std::complex<double> beta,
                                           /* */ std::complex<double>* c, KK_INT ldc) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  cblas_zherk(acc_order, acc_uplo, acc_trans, n, k, alpha.real(), a, lda, beta.real(), c, ldc);
#else
  F77_FUNC_ZHERK(&uplo, &transa, &n, &k, &alpha, (const std::complex<double>*)a, &lda, &beta, (std::complex<double>*)c,
                 &ldc);
#endif
}
template <>
void HostBlas<std::complex<double> >::trmm(const char side, const char uplo, const char transa, const char diag,
                                           KK_INT m, KK_INT n, const std::complex<double> alpha,
                                           const std::complex<double>* a, KK_INT lda,
                                           /* */ std::complex<double>* b, KK_INT ldb) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  CBLAS_SIDE acc_side        = kk_to_accelerate_side(side);
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_ztrmm(acc_order, acc_side, acc_uplo, acc_trans, acc_diag, m, n, &alpha, a, lda, b, ldb);
#else
  F77_FUNC_ZTRMM(&side, &uplo, &transa, &diag, &m, &n, &alpha, (const std::complex<double>*)a, &lda,
                 (std::complex<double>*)b, &ldb);
#endif
}
template <>
void HostBlas<std::complex<double> >::trsm(const char side, const char uplo, const char transa, const char diag,
                                           KK_INT m, KK_INT n, const std::complex<double> alpha,
                                           const std::complex<double>* a, KK_INT lda,
                                           /* */ std::complex<double>* b, KK_INT ldb) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  enum CBLAS_ORDER acc_order = CblasColMajor;
  CBLAS_SIDE acc_side        = kk_to_accelerate_side(side);
  enum CBLAS_UPLO acc_uplo   = kk_to_accelerate_uplo(uplo);
  auto acc_trans             = kk_to_accelerate_trans(transa);
  enum CBLAS_DIAG acc_diag   = kk_to_accelerate_diag(diag);
  cblas_ztrsm(acc_order, acc_side, acc_uplo, acc_trans, acc_diag, m, n, &alpha, a, lda, b, ldb);
#else
  F77_FUNC_ZTRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, (const std::complex<double>*)a, &lda,
                 (std::complex<double>*)b, &ldb);
#endif
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS
