/// \file KokkosBlas_Host_tpl.cpp
/// \brief BLAS wrapper for host tpls
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_config.h"
#include "KokkosBlas_Host_tpl.hpp"

#if defined( KOKKOSKERNELS_ENABLE_TPL_BLAS )

/// Fortran headers
extern "C" {

  ///
  /// scal
  ///
  void F77_BLAS_MANGLE(sscal,SSCAL)( const int* N, 
                                     const float* alpha,
                                     /* */ float* x, const int* x_inc);
  void F77_BLAS_MANGLE(dscal,DSCAL)( const int* N, 
                                     const double* alpha,
                                     /* */ double* x, const int* x_inc);
  void F77_BLAS_MANGLE(cscal,CSCAL)( const int* N, 
                                     const std::complex<float>* alpha,
                                     /* */ std::complex<float>* x, const int* x_inc);
  void F77_BLAS_MANGLE(zscal,ZSCAL)( const int* N, 
                                     const std::complex<double>* alpha,
                                     /* */ std::complex<double>* x, const int* x_inc);
  
  ///
  /// max
  ///
  int F77_BLAS_MANGLE(isamax,ISAMAX)( const int* N, const float* x, const int* x_inc);
  int F77_BLAS_MANGLE(idamax,IDAMAX)( const int* N, const double* x, const int* x_inc);
  int F77_BLAS_MANGLE(icamax,ICAMAX)( const int* N, const std::complex<float>* x, const int* x_inc);
  int F77_BLAS_MANGLE(izamax,IZAMAX)( const int* N, const std::complex<double>* x, const int* x_inc);

  
  ///
  /// nrm2
  ///
  float  F77_BLAS_MANGLE(snrm2, SNRM2 )( const int* N, const float* x, const int* x_inc);
  double F77_BLAS_MANGLE(dnrm2, DNRM2 )( const int* N, const double* x, const int* x_inc);
  float  F77_BLAS_MANGLE(scnrm2,SCNRM2)( const int* N, const std::complex<float>* x, const int* x_inc);
  double F77_BLAS_MANGLE(dznrm2,DZNRM2)( const int* N, const std::complex<double>* x, const int* x_inc);
  
  ///
  /// sum
  ///  
  float  F77_BLAS_MANGLE(sasum, SASUM )( const int* N, const float* x, const int* x_inc);
  double F77_BLAS_MANGLE(dasum, DASUM )( const int* N, const double* x, const int* x_inc);
  float  F77_BLAS_MANGLE(scasum,SCASUM)( const int* N, const std::complex<float>* x, const int* x_inc);
  double F77_BLAS_MANGLE(dzasum,DZASUM)( const int* N, const std::complex<double>* x, const int* x_inc);
    
  ///
  /// dot 
  ///
  float  F77_BLAS_MANGLE(sdot,SDOT)( const int* N, const float* x, const int* x_inc,
                                     const float* y, const int* y_inc);
  double F77_BLAS_MANGLE(ddot,DDOT)( const int* N, const double* x, const int* x_inc,
                                     const double* y, const int* y_inc);
# if defined( KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX )
  std::complex<float>  F77_BLAS_MANGLE(cdotu,CDOTU)( const int* N, 
                                                     const std::complex<float>* x, const int* x_inc,
                                                     const std::complex<float>* y, const int* y_inc);   
  std::complex<double> F77_BLAS_MANGLE(zdotu,ZDOTU)( const int* N, 
                                                     const std::complex<double>* x, const int* x_inc,
                                                     const std::complex<double>* y, const int* y_inc);
  std::complex<float>  F77_BLAS_MANGLE(cdotc,CDOTC)( const int* N, 
                                                     const std::complex<float>* x, const int* x_inc,
                                                     const std::complex<float>* y, const int* y_inc);  
  std::complex<double> F77_BLAS_MANGLE(zdotc,ZDOTC)( const int* N, 
                                                     const std::complex<double>* x, const int* x_inc,
                                                     const std::complex<double>* y, const int* y_inc);
# else 
  void F77_BLAS_MANGLE(cdotu,CDOTU)( std::complex<float> *res, 
                                     const int* N, 
                                     const std::complex<float>* x, const int* x_inc,
                                     const std::complex<float>* y, const int* y_inc);   
  void F77_BLAS_MANGLE(zdotu,ZDOTU)( std::complex<double> *res, 
                                     const int* N, 
                                     const std::complex<double>* x, const int* x_inc,
                                     const std::complex<double>* y, const int* y_inc);
  void F77_BLAS_MANGLE(cdotc,CDOTC)( std::complex<float> *res, 
                                     const int* N, 
                                     const std::complex<float>* x, const int* x_inc,
                                     const std::complex<float>* y, const int* y_inc); 
  void F77_BLAS_MANGLE(zdotc,ZDOTC)( std::complex<double> *res, 
                                     const int* N, 
                                     const std::complex<double>* x, const int* x_inc,
                                     const std::complex<double>* y, const int* y_inc);
# endif
  
  ///
  /// axpy
  ///
  void F77_BLAS_MANGLE(saxpy,SAXPY)( const int* N, 
                                     const float* alpha,
                                     const float* x, const int* x_inc,
                                     /* */ float* y, const int* y_inc);
  void F77_BLAS_MANGLE(daxpy,DAXPY)( const int* N, 
                                     const double* alpha,
                                     const double* x, const int* x_inc,
                                     /* */ double* y, const int* y_inc);
  void F77_BLAS_MANGLE(caxpy,CAXPY)( const int* N, 
                                     const std::complex<float>* alpha,
                                     const std::complex<float>* x, const int* x_inc,
                                     /* */ std::complex<float>* y, const int* y_inc);
  void F77_BLAS_MANGLE(zaxpy,ZAXPY)( const int* N, 
                                     const std::complex<double>* alpha,
                                     const std::complex<double>* x, const int* x_inc,
                                     /* */ std::complex<double>* y, const int* y_inc);
  
  ///
  /// Gemv
  ///
  void F77_BLAS_MANGLE(sgemv,SGEMV)( const char*, 
                                     int*, int*, 
                                     const float*,
                                     const float*, int*,
                                     const float*, int*,
                                     const float*,
                                     /* */ float*, int* );
  void F77_BLAS_MANGLE(dgemv,DGEMV)( const char*,
                                     int*, int*, 
                                     const double*,
                                     const double*, int*,
                                     const double*, int*,
                                     const double*,
                                     /* */ double*, int* );
  void F77_BLAS_MANGLE(cgemv,CGEMV)( const char*,
                                     int*, int*, 
                                     const std::complex<float>*,
                                     const std::complex<float>*, int*,
                                     const std::complex<float>*, int*,
                                     const std::complex<float>*,
                                     /* */ std::complex<float>*, int* );
  void F77_BLAS_MANGLE(zgemv,ZGEMV)( const char*,
                                     int*, int*, 
                                     const std::complex<double>*,
                                     const std::complex<double>*, int*,
                                     const std::complex<double>*, int*,
                                     const std::complex<double>*,
                                     /* */ std::complex<double>*, int* );

  ///
  /// Trsv
  ///

  void F77_BLAS_MANGLE(strsv,STRSV)( const char*, const char*, const char*, 
                                     int*, 
                                     const float*, int*,
                                     /* */ float*, int* );
  void F77_BLAS_MANGLE(dtrsv,DTRSV)( const char*, const char*, const char*, 
                                     int*, 
                                     const double*, int*,
                                     /* */ double*, int* );
  void F77_BLAS_MANGLE(ctrsv,CTRSV)( const char*, const char*, const char*, 
                                     int*, 
                                     const std::complex<float>*, int*,
                                     /* */ std::complex<float>*, int* );
  void F77_BLAS_MANGLE(ztrsv,ZTRSV)( const char*, const char*, const char*, 
                                     int*, 
                                     const std::complex<double>*, int*,
                                     /* */ std::complex<double>*, int* );

  ///
  /// Gemm
  ///

  void F77_BLAS_MANGLE(sgemm,SGEMM)( const char*, const char*,
                                     int*, int*, int*,
                                     const float*,
                                     const float*, int*,
                                     const float*, int*,
                                     const float*,
                                     /* */ float*, int* );
  void F77_BLAS_MANGLE(dgemm,DGEMM)( const char*, const char*,
                                     int*, int*, int*,
                                     const double*,
                                     const double*, int*,
                                     const double*, int*,
                                     const double*,
                                     /* */ double*, int* );
  void F77_BLAS_MANGLE(cgemm,CGEMM)( const char*, const char*,
                                     int*, int*, int*,
                                     const std::complex<float>*,
                                     const std::complex<float>*, int*,
                                     const std::complex<float>*, int*,
                                     const std::complex<float>*,
                                     /* */ std::complex<float>*, int* );
  void F77_BLAS_MANGLE(zgemm,ZGEMM)( const char*, const char*,
                                     int*, int*, int*,
                                     const std::complex<double>*,
                                     const std::complex<double>*, int*,
                                     const std::complex<double>*, int*,
                                     const std::complex<double>*,
                                     /* */ std::complex<double>*, int* );

  ///
  /// Herk
  ///

  void F77_BLAS_MANGLE(ssyrk,SSYRK)( const char*, const char*,
                                     int*, int*, 
                                     const float*,
                                     const float*, int*,
                                     const float*,
                                     /* */ float*, int* );
  void F77_BLAS_MANGLE(dsyrk,DSYRK)( const char*, const char*,
                                     int*, int*, 
                                     const double*,
                                     const double*, int*,
                                     const double*,
                                     /* */ double*, int* );
  void F77_BLAS_MANGLE(cherk,CHERK)( const char*, const char*,
                                     int*, int*, 
                                     const std::complex<float>*,
                                     const std::complex<float>*, int*,
                                     const std::complex<float>*,
                                     /* */ std::complex<float>*, int* );
  void F77_BLAS_MANGLE(zherk,ZHERK)( const char*, const char*,
                                     int*, int*, 
                                     const std::complex<double>*,
                                     const std::complex<double>*, int*,
                                     const std::complex<double>*,
                                     /* */ std::complex<double>*, int* );

  ///
  /// Trmm
  ///

  void F77_BLAS_MANGLE(strmm,STRMM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const float*,
                                     const float*, int*,
                                     /* */ float*, int* );
  void F77_BLAS_MANGLE(dtrmm,DTRMM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const double*,
                                     const double*, int*,
                                     /* */ double*, int* );
  void F77_BLAS_MANGLE(ctrmm,CTRMM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const std::complex<float>*,
                                     const std::complex<float>*, int*,
                                     /* */ std::complex<float>*, int* );
  void F77_BLAS_MANGLE(ztrmm,ZTRMM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const std::complex<double>*,
                                     const std::complex<double>*, int*,
                                     /* */ std::complex<double>*, int* );

  ///
  /// Trsm
  ///

  void F77_BLAS_MANGLE(strsm,STRSM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const float*,
                                     const float*, int*,
                                     /* */ float*, int* );
  void F77_BLAS_MANGLE(dtrsm,DTRSM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const double*,
                                     const double*, int*,
                                     /* */ double*, int* );
  void F77_BLAS_MANGLE(ctrsm,CTRSM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const std::complex<float>*,
                                     const std::complex<float>*, int*,
                                     /* */ std::complex<float>*, int* );
  void F77_BLAS_MANGLE(ztrsm,ZTRSM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const std::complex<double>*,
                                     const std::complex<double>*, int*,
                                     /* */ std::complex<double>*, int* );

  ///
  /// Gesv
  ///

  void F77_BLAS_MANGLE(sgesv,SGESV)( int*, int*,
                                     float*, int*, int*,
                                     float*, int*, int* );
  void F77_BLAS_MANGLE(dgesv,DGESV)( int*, int*,
                                     double*, int*, int*,
                                     double*, int*, int* );
  void F77_BLAS_MANGLE(cgesv,CGESV)( int*, int*,
                                     std::complex<float>*, int*, int*,
                                     std::complex<float>*, int*, int* );
  void F77_BLAS_MANGLE(zgesv,ZGESV)( int*, int*,
                                     std::complex<double>*, int*, int*,
                                     std::complex<double>*, int*, int* );
}



  void F77_BLAS_MANGLE(sscal,SSCAL)( const int* N, 
                                     const float* alpha,
                                     /* */ float* x, const int* x_inc);
  void F77_BLAS_MANGLE(dscal,DSCAL)( const int* N, 
                                     const double* alpha,
                                     /* */ double* x, const int* x_inc);
  void F77_BLAS_MANGLE(cscal,CSCAL)( const int* N, 
                                     const std::complex<float>* alpha,
                                     /* */ std::complex<float>* x, const int* x_inc);
  void F77_BLAS_MANGLE(zscal,ZSCAL)( const int* N, 
                                     const std::complex<double>* alpha,
                                     /* */ std::complex<double>* x, const int* x_inc);
  
#define F77_FUNC_SSCAL F77_BLAS_MANGLE(sscal,SSCAL)
#define F77_FUNC_DSCAL F77_BLAS_MANGLE(dscal,DSCAL)
#define F77_FUNC_CSCAL F77_BLAS_MANGLE(cscal,CSCAL)
#define F77_FUNC_ZSCAL F77_BLAS_MANGLE(zscal,ZSCAL)

#define F77_FUNC_ISAMAX F77_BLAS_MANGLE(isamax,ISAMAX)
#define F77_FUNC_IDAMAX F77_BLAS_MANGLE(idamax,IDAMAX)
#define F77_FUNC_ICAMAX F77_BLAS_MANGLE(icamax,ICAMAX)
#define F77_FUNC_IZAMAX F77_BLAS_MANGLE(izamax,IZAMAX)
  
#define F77_FUNC_SNRM2  F77_BLAS_MANGLE(snrm2,  SNRM2 )
#define F77_FUNC_DNRM2  F77_BLAS_MANGLE(dnrm2,  DNRM2 )
#define F77_FUNC_SCNRM2 F77_BLAS_MANGLE(scnrm2, SCNRM2)
#define F77_FUNC_DZNRM2 F77_BLAS_MANGLE(dznrm2, DZNRM2)
      
#define F77_FUNC_SASUM  F77_BLAS_MANGLE(sasum,  SASUM )
#define F77_FUNC_DASUM  F77_BLAS_MANGLE(dasum,  DASUM )
#define F77_FUNC_SCASUM F77_BLAS_MANGLE(scasum, SCASUM)
#define F77_FUNC_DZASUM F77_BLAS_MANGLE(dzasum, DZASUM)

#define F77_FUNC_SDOT  F77_BLAS_MANGLE(sdot,SDOT)
#define F77_FUNC_DDOT  F77_BLAS_MANGLE(ddot,DDOT)
#define F77_FUNC_CDOTU F77_BLAS_MANGLE(cdotu,CDOTU)
#define F77_FUNC_ZDOTU F77_BLAS_MANGLE(zdotu,ZDOTU)
#define F77_FUNC_CDOTC F77_BLAS_MANGLE(cdotc,CDOTC)
#define F77_FUNC_ZDOTC F77_BLAS_MANGLE(zdotc,ZDOTC)

#define F77_FUNC_SAXPY F77_BLAS_MANGLE(saxpy,SAXPY)
#define F77_FUNC_DAXPY F77_BLAS_MANGLE(daxpy,DAXPY)
#define F77_FUNC_CAXPY F77_BLAS_MANGLE(caxpy,CAXPY)
#define F77_FUNC_ZAXPY F77_BLAS_MANGLE(zaxpy,ZAXPY)

#define F77_FUNC_SGEMV F77_BLAS_MANGLE(sgemv,SGEMV)
#define F77_FUNC_DGEMV F77_BLAS_MANGLE(dgemv,DGEMV)
#define F77_FUNC_CGEMV F77_BLAS_MANGLE(cgemv,CGEMV)
#define F77_FUNC_ZGEMV F77_BLAS_MANGLE(zgemv,ZGEMV)

#define F77_FUNC_STRSV F77_BLAS_MANGLE(strsv,STRSV)
#define F77_FUNC_DTRSV F77_BLAS_MANGLE(dtrsv,DTRSV)
#define F77_FUNC_CTRSV F77_BLAS_MANGLE(ctrsv,CTRSV)
#define F77_FUNC_ZTRSV F77_BLAS_MANGLE(ztrsv,ZTRSV)

#define F77_FUNC_SGEMM F77_BLAS_MANGLE(sgemm,SGEMM)
#define F77_FUNC_DGEMM F77_BLAS_MANGLE(dgemm,DGEMM)
#define F77_FUNC_CGEMM F77_BLAS_MANGLE(cgemm,CGEMM)
#define F77_FUNC_ZGEMM F77_BLAS_MANGLE(zgemm,ZGEMM)

#define F77_FUNC_SSYRK F77_BLAS_MANGLE(ssyrk,SSYRK)
#define F77_FUNC_DSYRK F77_BLAS_MANGLE(dsyrk,DSYRK)
#define F77_FUNC_CHERK F77_BLAS_MANGLE(cherk,CHERK)
#define F77_FUNC_ZHERK F77_BLAS_MANGLE(zherk,ZHERK)

#define F77_FUNC_STRMM F77_BLAS_MANGLE(strmm,STRMM)
#define F77_FUNC_DTRMM F77_BLAS_MANGLE(dtrmm,DTRMM)
#define F77_FUNC_CTRMM F77_BLAS_MANGLE(ctrmm,CTRMM)
#define F77_FUNC_ZTRMM F77_BLAS_MANGLE(ztrmm,ZTRMM)

#define F77_FUNC_STRSM F77_BLAS_MANGLE(strsm,STRSM)
#define F77_FUNC_DTRSM F77_BLAS_MANGLE(dtrsm,DTRSM)
#define F77_FUNC_CTRSM F77_BLAS_MANGLE(ctrsm,CTRSM)
#define F77_FUNC_ZTRSM F77_BLAS_MANGLE(ztrsm,ZTRSM)

#define F77_FUNC_SGESV F77_BLAS_MANGLE(sgesv,SGESV)
#define F77_FUNC_DGESV F77_BLAS_MANGLE(dgesv,DGESV)
#define F77_FUNC_CGESV F77_BLAS_MANGLE(cgesv,CGESV)
#define F77_FUNC_ZGESV F77_BLAS_MANGLE(zgesv,ZGESV)

namespace KokkosBlas {
  namespace Impl {

    ///
    /// float
    ///

    template<>
    void
    HostBlas<float>::scal(int n, 
                             const float alpha,
                             /* */ float *x, int x_inc) {
      F77_FUNC_SSCAL(&n, &alpha, x, &x_inc);
    }
    template<>
    int
    HostBlas<float>::iamax(int n, 
                              const float *x, int x_inc) {
      return F77_FUNC_ISAMAX(&n, x, &x_inc);
    }
    template<>
    float
    HostBlas<float>::nrm2(int n, 
                             const float *x, int x_inc) {
      return F77_FUNC_SNRM2(&n, x, &x_inc);
    }
    template<>
    float
    HostBlas<float>::asum(int n, 
                             const float *x, int x_inc) {
      return F77_FUNC_SASUM(&n, x, &x_inc);
    }
    template<>
    float
    HostBlas<float>::dot(int n, 
                            const float *x, int x_inc,
                            const float *y, int y_inc) {
      return F77_FUNC_SDOT(&n, x, &x_inc, y, &y_inc);
    }
    template<>
    void
    HostBlas<float>::axpy(int n,
                             const float alpha,
                             const float *x, int x_inc,
                             /* */ float *y, int y_inc) {
      F77_FUNC_SAXPY(&n, &alpha, x, &x_inc, y, &y_inc);
    }   
    template<>
    void
    HostBlas<float>::gemv(const char trans, 
                      int m, int n, 
                      const float alpha, 
                      const float *a, int lda,
                      const float *b, int ldb,
                      const float beta,
                      /* */ float *c, int ldc) {
      F77_FUNC_SGEMV(&trans, 
                     &m, &n,
                     &alpha,
                     a, &lda,
                     b, &ldb,
                     &beta,
                     c, &ldc);
    }
    template<>
    void
    HostBlas<float>::trsv(const char uplo, const char transa, const char diag, 
                      int m, 
                      const float *a, int lda,
                      /* */ float *b, int ldb) {
      F77_FUNC_STRSV(&uplo, &transa, &diag,
                     &m,
                     a, &lda,
                     b, &ldb);
    }
    template<>
    void 
    HostBlas<float>::gemm(const char transa, const char transb, 
                      int m, int n, int k,
                      const float alpha, 
                      const float *a, int lda,
                      const float *b, int ldb,
                      const float beta,
                      /* */ float *c, int ldc) {
      F77_FUNC_SGEMM(&transa, &transb,
                     &m, &n, &k,
                     &alpha,
                     a, &lda,
                     b, &ldb,
                     &beta,
                     c, &ldc);
    }
    template<>
    void 
    HostBlas<float>::herk(const char transa, const char transb, 
                      int n, int k,
                      const float alpha, 
                      const float *a, int lda,
                      const float beta,
                      /* */ float *c, int ldc) {
      F77_FUNC_SSYRK(&transa, &transb,
                     &n, &k,
                     &alpha,
                     a, &lda,
                     &beta,
                     c, &ldc);
    }
    template<>
    void 
    HostBlas<float>::trmm(const char side, const char uplo, const char transa, const char diag,
                      int m, int n, 
                      const float alpha, 
                      const float *a, int lda,
                      /* */ float *b, int ldb) {
      F77_FUNC_STRMM(&side, &uplo, &transa, &diag,
                     &m, &n,
                     &alpha,
                     a, &lda,
                     b, &ldb);
    }
    template<>
    void 
    HostBlas<float>::trsm(const char side, const char uplo, const char transa, const char diag,
                      int m, int n, 
                      const float alpha, 
                      const float *a, int lda,
                      /* */ float *b, int ldb) {
      F77_FUNC_STRSM(&side, &uplo, &transa, &diag,
                     &m, &n,
                     &alpha,
                     a, &lda,
                     b, &ldb);
    }
    template<>
    void 
    HostBlas<float>::gesv(int n, int rhs,
                          float *a, int lda, int *ipiv,
                          float *b, int ldb, int info) {
      F77_FUNC_SGESV(&n, &rhs,
                     a, &lda, ipiv,
                     b, &ldb, &info);
    }

    ///
    /// double
    ///

    template<>
    void
    HostBlas<double>::scal(int n, 
                             const double alpha,
                             /* */ double *x, int x_inc) {
      F77_FUNC_DSCAL(&n, &alpha, x, &x_inc);
    }
    template<>
    int
    HostBlas<double>::iamax(int n, 
                               const double *x, int x_inc) {
      return F77_FUNC_IDAMAX(&n, x, &x_inc);
    }
    template<>
    double
    HostBlas<double>::nrm2(int n, 
                              const double *x, int x_inc) {
      return F77_FUNC_DNRM2(&n, x, &x_inc);
    }
    template<>
    double
    HostBlas<double>::asum(int n, 
                              const double *x, int x_inc) {
      return F77_FUNC_DASUM(&n, x, &x_inc);
    }
    template<>
    double
    HostBlas<double>::dot(int n, 
                            const double *x, int x_inc,
                            const double *y, int y_inc) {
      return F77_FUNC_DDOT(&n, x, &x_inc, y, &y_inc);
    }
    template<>
    void
    HostBlas<double>::axpy(int n,
                             const double alpha,
                             const double *x, int x_inc,
                             /* */ double *y, int y_inc) {
      F77_FUNC_DAXPY(&n, &alpha, x, &x_inc, y, &y_inc);
    }     
    template<>
    void 
    HostBlas<double>::gemv(const char trans, 
                       int m, int n, 
                       const double alpha, 
                       const double *a, int lda,
                       const double *b, int ldb,
                       const double beta,
                       /* */ double *c, int ldc) {
      F77_FUNC_DGEMV(&trans, 
                     &m, &n,
                     &alpha,
                     a, &lda,
                     b, &ldb,
                     &beta,
                     c, &ldc);
    }
    template<>
    void 
    HostBlas<double>::trsv(const char uplo, const char transa, const char diag, 
                       int m, 
                       const double *a, int lda,
                       /* */ double *b, int ldb) {
      F77_FUNC_DTRSV(&uplo, &transa, &diag,
                     &m,
                     a, &lda,
                     b, &ldb);
    }
    template<>
    void 
    HostBlas<double>::gemm(const char transa, const char transb, 
                       int m, int n, int k,
                       const double alpha, 
                       const double *a, int lda,
                       const double *b, int ldb,
                       const double beta,
                       /* */ double *c, int ldc) {
      F77_FUNC_DGEMM(&transa, &transb,
                     &m, &n, &k,
                     &alpha,
                     a, &lda,
                     b, &ldb,
                     &beta,
                     c, &ldc);
    }
    template<>
    void 
    HostBlas<double>::herk(const char transa, const char transb, 
                       int n, int k,
                       const double alpha, 
                       const double *a, int lda,
                       const double beta,
                       /* */ double *c, int ldc) {
      F77_FUNC_DSYRK(&transa, &transb,
                     &n, &k,
                     &alpha,
                     a, &lda,
                     &beta,
                     c, &ldc);
    }
    template<>
    void 
    HostBlas<double>::trmm(const char side, const char uplo, const char transa, const char diag,
                       int m, int n, 
                       const double alpha, 
                       const double *a, int lda,
                       /* */ double *b, int ldb) {
      F77_FUNC_DTRMM(&side, &uplo, &transa, &diag,
                     &m, &n,
                     &alpha,
                     a, &lda,
                     b, &ldb);
    }
    template<>
    void 
    HostBlas<double>::trsm(const char side, const char uplo, const char transa, const char diag,
                       int m, int n, 
                       const double alpha, 
                       const double *a, int lda,
                       /* */ double *b, int ldb) {
      F77_FUNC_DTRSM(&side, &uplo, &transa, &diag,
                     &m, &n,
                     &alpha,
                     a, &lda,
                     b, &ldb);
    }
    template<>
    void 
    HostBlas<double>::gesv(int n, int rhs,
                          double *a, int lda, int *ipiv,
                          double *b, int ldb, int info) {
      F77_FUNC_DGESV(&n, &rhs,
                     a, &lda, ipiv,
                     b, &ldb, &info);
    }

    /// 
    /// std::complex<float>
    ///

    template<>
    void
    HostBlas<std::complex<float> >::scal(int n, 
                                            const std::complex<float> alpha,
                                            /* */ std::complex<float> *x, int x_inc) {
      F77_FUNC_CSCAL(&n, &alpha, x, &x_inc);
    }
    template<>
    int
    HostBlas<std::complex<float> >::iamax(int n, 
                                             const std::complex<float> *x, int x_inc) {
      return F77_FUNC_ICAMAX(&n, x, &x_inc);
    }
    template<>
    float
    HostBlas<std::complex<float> >::nrm2(int n, 
                                            const std::complex<float> *x, int x_inc) {
      return F77_FUNC_SCNRM2(&n, x, &x_inc);
    }
    template<>
    float
    HostBlas<std::complex<float> >::asum(int n, 
                                            const std::complex<float> *x, int x_inc) {
      return F77_FUNC_SCASUM(&n, x, &x_inc);
    }
    template<>
    std::complex<float>
    HostBlas<std::complex<float> >::dot(int n, 
                                           const std::complex<float> *x, int x_inc,
                                           const std::complex<float> *y, int y_inc) {
# if defined( KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX )
      return F77_FUNC_CDOTC(&n, x, &x_inc, y, &y_inc);
# else
      std::complex<float> res;
      F77_FUNC_CDOTC(&res, &n, x, &x_inc, y, &y_inc);
      return res;
# endif
    }
    template<>
    void
    HostBlas<std::complex<float> >::axpy(int n,
                                            const std::complex<float> alpha,
                                            const std::complex<float> *x, int x_inc,
                                            /* */ std::complex<float> *y, int y_inc) {
      F77_FUNC_CAXPY(&n, &alpha, x, &x_inc, y, &y_inc);
    }     
    
    template<>
    void 
    HostBlas<std::complex<float> >::gemv(const char trans, 
                                        int m, int n, 
                                        const std::complex<float> alpha, 
                                        const std::complex<float> *a, int lda,
                                        const std::complex<float> *b, int ldb,
                                        const std::complex<float> beta,
                                        /* */ std::complex<float> *c, int ldc) {
      F77_FUNC_CGEMV(&trans, 
                     &m, &n,
                     &alpha,
                     (const std::complex<float>*)a, &lda,
                     (const std::complex<float>*)b, &ldb,
                     &beta,
                     (      std::complex<float>*)c, &ldc);
    }
    template<>
    void 
    HostBlas<std::complex<float> >::trsv(const char uplo, const char transa, const char diag, 
                                        int m, 
                                        const std::complex<float> *a, int lda,
                                        /* */ std::complex<float> *b, int ldb) {
      F77_FUNC_CTRSV(&uplo, &transa, &diag,
                     &m,
                     (const std::complex<float>*)a, &lda,
                     (      std::complex<float>*)b, &ldb);
    }
    template<>
    void 
    HostBlas<std::complex<float> >::gemm(const char transa, const char transb, 
                                        int m, int n, int k,
                                        const std::complex<float> alpha, 
                                        const std::complex<float> *a, int lda,
                                        const std::complex<float> *b, int ldb,
                                        const std::complex<float> beta,
                                        /* */ std::complex<float> *c, int ldc) {
      F77_FUNC_CGEMM(&transa, &transb,
                     &m, &n, &k,
                     &alpha,
                     (const std::complex<float>*)a, &lda,
                     (const std::complex<float>*)b, &ldb,
                     &beta,
                     (      std::complex<float>*)c, &ldc);
    }
    template<>
    void 
    HostBlas<std::complex<float> >::herk(const char transa, const char transb, 
                                        int n, int k,
                                        const std::complex<float> alpha, 
                                        const std::complex<float> *a, int lda,
                                        const std::complex<float> beta,
                                        /* */ std::complex<float> *c, int ldc) {
      F77_FUNC_CHERK(&transa, &transb,
                     &n, &k,
                     &alpha,
                     (const std::complex<float>*)a, &lda,
                     &beta,
                     (      std::complex<float>*)c, &ldc);
    }
    template<>
    void 
    HostBlas<std::complex<float> >::trmm(const char side, const char uplo, const char transa, const char diag,
                                        int m, int n, 
                                        const std::complex<float> alpha, 
                                        const std::complex<float> *a, int lda,
                                        /* */ std::complex<float> *b, int ldb) {
      F77_FUNC_CTRMM(&side, &uplo, &transa, &diag,
                     &m, &n,
                     &alpha,
                     (const std::complex<float>*)a, &lda,
                     (      std::complex<float>*)b, &ldb);
    }
    template<>
    void 
    HostBlas<std::complex<float> >::trsm(const char side, const char uplo, const char transa, const char diag,
                                        int m, int n, 
                                        const std::complex<float> alpha, 
                                        const std::complex<float> *a, int lda,
                                        /* */ std::complex<float> *b, int ldb) {
      F77_FUNC_CTRSM(&side, &uplo, &transa, &diag,
                     &m, &n,
                     &alpha,
                     (const std::complex<float>*)a, &lda,
                     (      std::complex<float>*)b, &ldb);
    }
    template<>
    void 
    HostBlas<std::complex<float> >::gesv(int n, int rhs,
                                         std::complex<float> *a, int lda, int *ipiv,
                                         std::complex<float> *b, int ldb, int info) {
      F77_FUNC_CGESV(&n, &rhs,
                     a, &lda, ipiv,
                     b, &ldb, &info);
    }
    
    ///
    /// std::complex<double>
    ///
    

    template<>
    void
    HostBlas<std::complex<double> >::scal(int n, 
                                            const std::complex<double> alpha,
                                            /* */ std::complex<double> *x, int x_inc) {
      F77_FUNC_ZSCAL(&n, &alpha, x, &x_inc);
    }
    template<>
    int
    HostBlas<std::complex<double> >::iamax(int n, 
                                              const std::complex<double> *x, int x_inc) {
      return F77_FUNC_IZAMAX(&n, x, &x_inc);
    }
    template<>
    double
    HostBlas<std::complex<double> >::nrm2(int n, 
                                            const std::complex<double> *x, int x_inc) {
      return F77_FUNC_DZNRM2(&n, x, &x_inc);
    }
    template<>
    double
    HostBlas<std::complex<double> >::asum(int n, 
                                            const std::complex<double> *x, int x_inc) {
      return F77_FUNC_DZASUM(&n, x, &x_inc);
    }
    template<>
    std::complex<double>
    HostBlas<std::complex<double> >::dot(int n, 
                                           const std::complex<double> *x, int x_inc,
                                           const std::complex<double> *y, int y_inc) {
# if defined( KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX )
      return F77_FUNC_ZDOTC(&n, x, &x_inc, y, &y_inc);
# else
      std::complex<double> res;
      F77_FUNC_ZDOTC(&res, &n, x, &x_inc, y, &y_inc);
      return res;
# endif
    }
    template<>
    void
    HostBlas<std::complex<double> >::axpy(int n,
                                            const std::complex<double> alpha,
                                            const std::complex<double> *x, int x_inc,
                                            /* */ std::complex<double> *y, int y_inc) {
      F77_FUNC_ZAXPY(&n, &alpha, x, &x_inc, y, &y_inc);
    }     

    template<>
    void 
    HostBlas<std::complex<double> >::gemv(const char trans, 
                                         int m, int n, 
                                         const std::complex<double> alpha, 
                                         const std::complex<double> *a, int lda,
                                         const std::complex<double> *b, int ldb,
                                         const std::complex<double> beta,
                                         /* */ std::complex<double> *c, int ldc) {
      F77_FUNC_ZGEMV(&trans, 
                     &m, &n,
                     &alpha,
                     (const std::complex<double>*)a, &lda,
                     (const std::complex<double>*)b, &ldb,
                     &beta,
                     (      std::complex<double>*)c, &ldc);
    }
    template<>
    void 
    HostBlas<std::complex<double> >::trsv(const char uplo, const char transa, const char diag, 
                                         int m, 
                                         const std::complex<double> *a, int lda,
                                         /* */ std::complex<double> *b, int ldb) {
      F77_FUNC_ZTRSV(&uplo, &transa, &diag,
                     &m,
                     (const std::complex<double>*)a, &lda,
                     (      std::complex<double>*)b, &ldb);
    }
    template<>
    void 
    HostBlas<std::complex<double> >::gemm(const char transa, const char transb, 
                                         int m, int n, int k,
                                         const std::complex<double> alpha, 
                                         const std::complex<double> *a, int lda,
                                         const std::complex<double> *b, int ldb,
                                         const std::complex<double> beta,
                                         /* */ std::complex<double> *c, int ldc) {
      F77_FUNC_ZGEMM(&transa, &transb,
                     &m, &n, &k,
                     &alpha,
                     (const std::complex<double>*)a, &lda,
                     (const std::complex<double>*)b, &ldb,
                     &beta,
                     (      std::complex<double>*)c, &ldc);
    }
    template<>
    void 
    HostBlas<std::complex<double> >::herk(const char transa, const char transb, 
                                         int n, int k,
                                         const std::complex<double> alpha, 
                                         const std::complex<double> *a, int lda,
                                         const std::complex<double> beta,
                                         /* */ std::complex<double> *c, int ldc) {
      F77_FUNC_ZHERK(&transa, &transb,
                     &n, &k,
                     &alpha,
                     (const std::complex<double>*)a, &lda,
                     &beta,
                     (      std::complex<double>*)c, &ldc);
    }
    template<>
    void 
    HostBlas<std::complex<double> >::trmm(const char side, const char uplo, const char transa, const char diag,
                                         int m, int n, 
                                         const std::complex<double> alpha, 
                                         const std::complex<double> *a, int lda,
                                         /* */ std::complex<double> *b, int ldb) {
      F77_FUNC_ZTRMM(&side, &uplo, &transa, &diag,
                     &m, &n,
                     &alpha,
                     (const std::complex<double>*)a, &lda,
                     (      std::complex<double>*)b, &ldb);
    }
    template<>
    void 
    HostBlas<std::complex<double> >::trsm(const char side, const char uplo, const char transa, const char diag,
                                         int m, int n, 
                                         const std::complex<double> alpha, 
                                         const std::complex<double> *a, int lda,
                                         /* */ std::complex<double> *b, int ldb) {
      F77_FUNC_ZTRSM(&side, &uplo, &transa, &diag,
                     &m, &n,
                     &alpha,
                     (const std::complex<double>*)a, &lda,
                     (      std::complex<double>*)b, &ldb);
    }
    template<>
    void 
    HostBlas<std::complex<double> >::gesv(int n, int rhs,
                                         std::complex<double> *a, int lda, int *ipiv,
                                         std::complex<double> *b, int ldb, int info) {
      F77_FUNC_ZGESV(&n, &rhs,
                     a, &lda, ipiv,
                     b, &ldb, &info);
    }

  }
}
#endif
