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
}

namespace Tacho {

  namespace Experimental {

#define F77_FUNC_SPOTRF F77_BLAS_MANGLE(spotrf,SPOTRF)
#define F77_FUNC_DPOTRF F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define F77_FUNC_CPOTRF F77_BLAS_MANGLE(cpotrf,CPOTRF)
#define F77_FUNC_ZPOTRF F77_BLAS_MANGLE(zpotrf,ZPOTRF)

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
    };

    template<typename T>
    struct Lapack<T, typename std::enable_if<std::is_same<T, std   ::complex<float> >::value ||
                                             std::is_same<T, Kokkos::complex<float> >::value, T>::type> {
      static inline
      void potrf(const char uplo,
                 int m,
                 T *a, int lda,
                 int *info) {
        F77_FUNC_CPOTRF(&uplo,
                        &m,
                        a, &lda,
                        info);
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
        F77_FUNC_ZPOTRF(&uplo,
                        &m,
                        a, &lda,
                        info);
      }
    };
  }
}

#endif
