#ifndef __TACHOEXP_BLAS_EXTERNAL_HPP__
#define __TACHOEXP_BLAS_EXTERNAL_HPP__

/// \file  TachoExp_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "ShyLU_NodeTacho_config.h"
#include "TachoExp_Util.hpp"


extern "C" {
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
                                     const Kokkos::complex<float>*,
                                     const Kokkos::complex<float>*, int*,
                                     const Kokkos::complex<float>*, int*,
                                     const Kokkos::complex<float>*,
                                     /* */ Kokkos::complex<float>*, int* );
  void F77_BLAS_MANGLE(zgemv,ZGEMV)( const char*,
                                     int*, int*, 
                                     const Kokkos::complex<double>*,
                                     const Kokkos::complex<double>*, int*,
                                     const Kokkos::complex<double>*, int*,
                                     const Kokkos::complex<double>*,
                                     /* */ Kokkos::complex<double>*, int* );

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
                                     const Kokkos::complex<float>*, int*,
                                     /* */ Kokkos::complex<float>*, int* );
  void F77_BLAS_MANGLE(ztrsv,ZTRSV)( const char*, const char*, const char*, 
                                     int*, 
                                     const Kokkos::complex<double>*, int*,
                                     /* */ Kokkos::complex<double>*, int* );

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
                                     const Kokkos::complex<float>*,
                                     const Kokkos::complex<float>*, int*,
                                     const Kokkos::complex<float>*, int*,
                                     const Kokkos::complex<float>*,
                                     /* */ Kokkos::complex<float>*, int* );
  void F77_BLAS_MANGLE(zgemm,ZGEMM)( const char*, const char*,
                                     int*, int*, int*,
                                     const Kokkos::complex<double>*,
                                     const Kokkos::complex<double>*, int*,
                                     const Kokkos::complex<double>*, int*,
                                     const Kokkos::complex<double>*,
                                     /* */ Kokkos::complex<double>*, int* );

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
                                     const Kokkos::complex<float>*,
                                     const Kokkos::complex<float>*, int*,
                                     const Kokkos::complex<float>*,
                                     /* */ Kokkos::complex<float>*, int* );
  void F77_BLAS_MANGLE(zherk,ZHERK)( const char*, const char*,
                                     int*, int*, 
                                     const Kokkos::complex<double>*,
                                     const Kokkos::complex<double>*, int*,
                                     const Kokkos::complex<double>*,
                                     /* */ Kokkos::complex<double>*, int* );

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
                                     const Kokkos::complex<float>*,
                                     const Kokkos::complex<float>*, int*,
                                     /* */ Kokkos::complex<float>*, int* );
  void F77_BLAS_MANGLE(ztrsm,ZTRSM)( const char*, const char*, const char*, const char*,
                                     int*, int*,
                                     const Kokkos::complex<double>*,
                                     const Kokkos::complex<double>*, int*,
                                     /* */ Kokkos::complex<double>*, int* );
}

namespace Tacho {

  namespace Experimental {

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

#define F77_FUNC_STRSM F77_BLAS_MANGLE(strsm,STRSM)
#define F77_FUNC_DTRSM F77_BLAS_MANGLE(dtrsm,DTRSM)
#define F77_FUNC_CTRSM F77_BLAS_MANGLE(ctrsm,CTRSM)
#define F77_FUNC_ZTRSM F77_BLAS_MANGLE(ztrsm,ZTRSM)

    
    template<typename T, typename dummy = T> 
    struct Blas;

    template<typename T>
    struct Blas<T, typename std::enable_if<std::is_same<T, float>::value, T>::type> {
      static inline
      void gemv(const char trans, 
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_SGEMV(&trans, 
                       &m, &n,
                       &alpha,
                       a, &lda,
                       b, &ldb,
                       &beta,
                       c, &ldc);
      }
      static inline
      void trsv(const char uplo, const char transa, const char diag, 
                int m, 
                const T *a, int lda,
                /* */ T *b, int ldb) {
        F77_FUNC_STRSV(&uplo, &transa, &diag,
                       &m,
                       a, &lda,
                       b, &ldb);
      }
      static inline
      void gemm(const char transa, const char transb, 
                int m, int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_SGEMM(&transa, &transb,
                       &m, &n, &k,
                       &alpha,
                       a, &lda,
                       b, &ldb,
                       &beta,
                       c, &ldc);
      }
      static inline
      void herk(const char transa, const char transb, 
                int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_SSYRK(&transa, &transb,
                       &n, &k,
                       &alpha,
                       a, &lda,
                       &beta,
                       c, &ldc);
      }
      static inline
      void trsm(const char side, const char uplo, const char transa, const char diag,
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                /* */ T *b, int ldb) {
        F77_FUNC_STRSM(&side, &uplo, &transa, &diag,
                       &m, &n,
                       &alpha,
                       a, &lda,
                       b, &ldb);
      }
    };

    template<typename T>
    struct Blas<T, typename std::enable_if<std::is_same<T, double>::value, T>::type> {
      static inline
      void gemv(const char trans, 
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_DGEMV(&trans, 
                       &m, &n,
                       &alpha,
                       a, &lda,
                       b, &ldb,
                       &beta,
                       c, &ldc);
      }
      static inline
      void trsv(const char uplo, const char transa, const char diag, 
                int m, 
                const T *a, int lda,
                /* */ T *b, int ldb) {
        F77_FUNC_DTRSV(&uplo, &transa, &diag,
                       &m,
                       a, &lda,
                       b, &ldb);
      }
      static inline
      void gemm(const char transa, const char transb, 
                int m, int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_DGEMM(&transa, &transb,
                       &m, &n, &k,
                       &alpha,
                       a, &lda,
                       b, &ldb,
                       &beta,
                       c, &ldc);
      }
      static inline
      void herk(const char transa, const char transb, 
                int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_DSYRK(&transa, &transb,
                       &n, &k,
                       &alpha,
                       a, &lda,
                       &beta,
                       c, &ldc);
      }
      static inline
      void trsm(const char side, const char uplo, const char transa, const char diag,
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                /* */ T *b, int ldb) {
        F77_FUNC_DTRSM(&side, &uplo, &transa, &diag,
                       &m, &n,
                       &alpha,
                       a, &lda,
                       b, &ldb);
      }
    };

    template<typename T>
    struct Blas<T, typename std::enable_if<std::is_same<T, std   ::complex<float> >::value || 
                                           std::is_same<T, Kokkos::complex<float> >::value, T>::type> {
      static inline
      void gemv(const char trans, 
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_CGEMV(&trans, 
                       &m, &n,
                       &alpha,
                       (const Kokkos::complex<float>*)a, &lda,
                       (const Kokkos::complex<float>*)b, &ldb,
                       &beta,
                       (      Kokkos::complex<float>*)c, &ldc);
      }
      static inline
      void trsv(const char uplo, const char transa, const char diag, 
                int m, 
                const T *a, int lda,
                /* */ T *b, int ldb) {
        F77_FUNC_CTRSV(&uplo, &transa, &diag,
                       &m,
                       (const Kokkos::complex<float>*)a, &lda,
                       (      Kokkos::complex<float>*)b, &ldb);
      }
      static inline
      void gemm(const char transa, const char transb, 
                int m, int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_CGEMM(&transa, &transb,
                       &m, &n, &k,
                       &alpha,
                       (const Kokkos::complex<float>*)a, &lda,
                       (const Kokkos::complex<float>*)b, &ldb,
                       &beta,
                       (      Kokkos::complex<float>*)c, &ldc);
      }
      static inline
      void herk(const char transa, const char transb, 
                int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_CHERK(&transa, &transb,
                       &n, &k,
                       &alpha,
                       (const Kokkos::complex<float>*)a, &lda,
                       &beta,
                       (      Kokkos::complex<float>*)c, &ldc);
      }
      static inline
      void trsm(const char side, const char uplo, const char transa, const char diag,
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                /* */ T *b, int ldb) {
        F77_FUNC_CTRSM(&side, &uplo, &transa, &diag,
                       &m, &n,
                       &alpha,
                       (const Kokkos::complex<float>*)a, &lda,
                       (      Kokkos::complex<float>*)b, &ldb);
      }
    };

    template<typename T>
    struct Blas<T, typename std::enable_if<std::is_same<T, std   ::complex<double> >::value || 
                                           std::is_same<T, Kokkos::complex<double> >::value, T>::type> {
      static inline
      void gemv(const char trans, 
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_ZGEMV(&trans, 
                       &m, &n,
                       &alpha,
                       (const Kokkos::complex<double>*)a, &lda,
                       (const Kokkos::complex<double>*)b, &ldb,
                       &beta,
                       (      Kokkos::complex<double>*)c, &ldc);
      }
      static inline
      void trsv(const char uplo, const char transa, const char diag, 
                int m, 
                const T *a, int lda,
                /* */ T *b, int ldb) {
        F77_FUNC_ZTRSV(&uplo, &transa, &diag,
                       &m,
                       (const Kokkos::complex<double>*)a, &lda,
                       (      Kokkos::complex<double>*)b, &ldb);
      }
      static inline
      void gemm(const char transa, const char transb, 
                int m, int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_ZGEMM(&transa, &transb,
                       &m, &n, &k,
                       &alpha,
                       (const Kokkos::complex<double>*)a, &lda,
                       (const Kokkos::complex<double>*)b, &ldb,
                       &beta,
                       (      Kokkos::complex<double>*)c, &ldc);
      }
      static inline
      void herk(const char transa, const char transb, 
                int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T beta,
                /* */ T *c, int ldc) {
        F77_FUNC_ZHERK(&transa, &transb,
                       &n, &k,
                       &alpha,
                       (const Kokkos::complex<double>*)a, &lda,
                       &beta,
                       (      Kokkos::complex<double>*)c, &ldc);
      }
      static inline
      void trsm(const char side, const char uplo, const char transa, const char diag,
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                /* */ T *b, int ldb) {
        F77_FUNC_ZTRSM(&side, &uplo, &transa, &diag,
                       &m, &n,
                       &alpha,
                       (const Kokkos::complex<double>*)a, &lda,
                       (      Kokkos::complex<double>*)b, &ldb);
      }
    };
  }
}

#endif
