#ifndef __TACHOEXP_LAPACK_EXTERNAL_HPP__
#define __TACHOEXP_LAPACK_EXTERNAL_HPP__

/// \file  TachoExp_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "ShyLUTacho_config.h"
#include "TachoExp_Util.hpp"

extern "C" {
  ///
  /// POTRF
  ///
#if !defined( TACHO_HAVE_MKL )
  void F77_BLAS_MANGLE(spotrf,SPOTRF)( const char*, 
                                       int*, 
                                       float*, int*,
                                       int* );
  void F77_BLAS_MANGLE(dpotrf,DPOTRF)( const char*, 
                                       int*,
                                       double*, int*,
                                       int*);
  void F77_BLAS_MANGLE(cpotrf,CPOTRF)( const char*, 
                                       int*,
                                       Kokkos::complex<float>*, int*,
                                       int*);
  void F77_BLAS_MANGLE(zpotrf,ZPOTRF)( const char*, 
                                       int*, 
                                       Kokkos::complex<double>*, int*,
                                       int*);

  void F77_BLAS_MANGLE(ssytrf,SSYTRF)( const char*, 
                                       int*, 
                                       float*, int*,
                                       int* );
  void F77_BLAS_MANGLE(dsytrf,DSYTRF)( const char*, 
                                       int*,
                                       double*, int*,
                                       int*);
  void F77_BLAS_MANGLE(csytrf,CSYTRF)( const char*, 
                                       int*,
                                       Kokkos::complex<float>*, int*,
                                       int*);
  void F77_BLAS_MANGLE(zsytrf,ZSYTRF)( const char*, 
                                       int*, 
                                       Kokkos::complex<double>*, int*,
                                       int*);

#endif
}

namespace Tacho {

  namespace Experimental {

#define F77_FUNC_SPOTRF F77_BLAS_MANGLE(spotrf,SPOTRF)
#define F77_FUNC_DPOTRF F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define F77_FUNC_CPOTRF F77_BLAS_MANGLE(cpotrf,CPOTRF)
#define F77_FUNC_ZPOTRF F77_BLAS_MANGLE(zpotrf,ZPOTRF)

#define F77_FUNC_SSYTRF F77_BLAS_MANGLE(ssytrf,SSYTRF)
#define F77_FUNC_DSYTRF F77_BLAS_MANGLE(dsytrf,DSYTRF)
#define F77_FUNC_CSYTRF F77_BLAS_MANGLE(csytrf,CSYTRF)
#define F77_FUNC_ZSYTRF F77_BLAS_MANGLE(zsytrf,ZSYTRF)

    template<typename T, typename dummy = T>
    struct Lapack;

    template<typename T>
    struct Lapack<T, typename std::enable_if<std::is_same<T, float>::value, T>::type> {
      static inline
      void potrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *info) {
        F77_FUNC_SPOTRF(&uplo,
                        &m,
                        a, &lda,
                        info);
      }
      static inline
      void sytrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *ipiv,
                 int *info) {
        F77_FUNC_SSYTRF(&uplo,
                        &m,
                        a, &lda,
                        ipiv,
                        info);
      }
    };

    template<typename T>
    struct Lapack<T, typename std::enable_if<std::is_same<T, double>::value, T>::type> {
      static inline
      void potrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *info) {
        F77_FUNC_DPOTRF(&uplo,
                        &m,
                        a, &lda,
                        info);
      }
      static inline
      void sytrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *ipiv,
                 int *info) {
        F77_FUNC_DSYTRF(&uplo,
                        &m,
                        a, &lda,
                        ipiv,
                        info);
      }

    };

    template<typename T>
    struct Lapack<T, typename std::enable_if<std::is_same<T, std   ::complex<float> >::value ||
                                             std::is_same<T, Kokkos::complex<float> >::value, T>::type> {
      static inline
      void potrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *info) {
#if defined( TACHO_HAVE_MKL )
        F77_FUNC_CPOTRF(&uplo,
                        &m,
                        (MKL_Complex8 *)a, &lda,
                        info);
#else
        F77_FUNC_CPOTRF(&uplo,
                        &m,
                        a, &lda,
                        info);
#endif
      }
      static inline
      void sytrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *ipiv,
                 int *info) {
#if defined( TACHO_HAVE_MKL )
        F77_FUNC_CSYTRF(&uplo,
                        &m,
                        (MKL_Complex8 *)a, &lda,
                        ipiv,
                        info);
#else
        F77_FUNC_CSYTRF(&uplo,
                        &m,
                        a, &lda,
                        ipiv,
                        info);
#endif
      }
    };

    template<typename T>
    struct Lapack<T, typename std::enable_if<std::is_same<T, std   ::complex<double> >::value ||
                                             std::is_same<T, Kokkos::complex<double> >::value, T>::type> {
      static inline
      void potrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *info) {
#if defined( TACHO_HAVE_MKL )
        F77_FUNC_ZPOTRF(&uplo,
                        &m,
                        (MKL_Complex16 *)a, &lda,
                        info);
#else
        F77_FUNC_ZPOTRF(&uplo,
                        &m,
                        a, &lda,
                        info);
#endif
      }
      static inline
      void sytrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *ipiv,
                 int *info) {
#if defined( TACHO_HAVE_MKL )
        F77_FUNC_ZSYTRF(&uplo,
                        &m,
                        (MKL_Complex16 *)a, &lda,
                        ipiv,
                        info);
#else
        F77_FUNC_ZSYTRF(&uplo,
                        &m,
                        a, &lda,
                        ipiv,
                        info);
#endif
      }
    };
  }
}

#endif
