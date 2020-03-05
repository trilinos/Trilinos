/// \file  Tacho_Lapack_External.hpp
/// \brief Lapack wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "ShyLU_NodeTacho_config.h"
#include "Tacho_Lapack_External.hpp"

#include "Kokkos_Core.hpp"

extern "C" {
  ///
  /// POTRF
  ///
  void F77_BLAS_MANGLE(spotrf,SPOTRF)( const char*, 
                                       const int*, 
                                       float*, const int*,
                                       int* );
  void F77_BLAS_MANGLE(dpotrf,DPOTRF)( const char*, 
                                       const int*,
                                       double*, const int*,
                                       int*);
  void F77_BLAS_MANGLE(cpotrf,CPOTRF)( const char*, 
                                       const int*,
                                       Kokkos::complex<float>*, const int*,
                                       int*);
  void F77_BLAS_MANGLE(zpotrf,ZPOTRF)( const char*, 
                                       const int*, 
                                       Kokkos::complex<double>*, const int*,
                                       int*);

  void F77_BLAS_MANGLE(ssytrf,SSYTRF)( const char*, 
                                       const int*, 
                                       float*, const int*,
                                       int*,
                                       float*, int*,
                                       int*);
  void F77_BLAS_MANGLE(dsytrf,DSYTRF)( const char*, 
                                       const int*,
                                       double*, const int*,
                                       int*,
                                       double*, int*,
                                       int*);
  void F77_BLAS_MANGLE(csytrf,CSYTRF)( const char*, 
                                       const int*,
                                       Kokkos::complex<float>*, const int*,
                                       int*,
                                       Kokkos::complex<float>*, int *,
                                       int*);
  void F77_BLAS_MANGLE(zsytrf,ZSYTRF)( const char*, 
                                       const int*, 
                                       Kokkos::complex<double>*, const int*,
                                       int*,
                                       Kokkos::complex<double>*, int*,
                                       int*);
}

namespace Tacho {

#define F77_FUNC_SPOTRF F77_BLAS_MANGLE(spotrf,SPOTRF)
#define F77_FUNC_DPOTRF F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define F77_FUNC_CPOTRF F77_BLAS_MANGLE(cpotrf,CPOTRF)
#define F77_FUNC_ZPOTRF F77_BLAS_MANGLE(zpotrf,ZPOTRF)

#define F77_FUNC_SSYTRF F77_BLAS_MANGLE(ssytrf,SSYTRF)
#define F77_FUNC_DSYTRF F77_BLAS_MANGLE(dsytrf,DSYTRF)
#define F77_FUNC_CSYTRF F77_BLAS_MANGLE(csytrf,CSYTRF)
#define F77_FUNC_ZSYTRF F77_BLAS_MANGLE(zsytrf,ZSYTRF)

  template<>
  int 
  Lapack<float>::potrf(const char uplo,
                       const int m,
                       float *a, const int lda,
                       int *info) {
    F77_FUNC_SPOTRF(&uplo,
                    &m,
                    a, &lda,
                    info);
    return 0;
  }
#if defined(KOKKOS_ENABLE_CUDA)
  template<>
  int 
  Lapack<float>::potrf_buffersize(cusolverDnHandle_t handle,
                                  const cublasFillMode_t uplo,
                                  const int m,
                                  float *a, const int lda,
                                  int *lwork) {
    const int r_val = cusolverDnSpotrf_bufferSize(handle,
                                                  uplo,
                                                  m, 
                                                  a, lda,
                                                  lwork);
    return r_val;
  }

  template<>
  int 
  Lapack<float>::potrf(cusolverDnHandle_t handle,
                       const cublasFillMode_t uplo,
                       const int m,
                       float *a, const int lda,
                       float *w, const int lwork,
                       int *dev) {
    const int r_val = cusolverDnSpotrf(handle,
                                       uplo,
                                       m, 
                                       a, lda,
                                       w, lwork,
                                       dev);
    return r_val;
  }
#endif

  template<>
  int 
  Lapack<float>::sytrf(const char uplo,
                       const int m,
                       float *a, const int lda,
                       int *ipiv,
                       float *work, int lwork,
                       int *info) {
    F77_FUNC_SSYTRF(&uplo,
                    &m,
                    a, &lda,
                    ipiv,
                    work, &lwork,
                    info);
    return 0;
  }
    
  template<>
  int
  Lapack<double>::potrf(const char uplo,
                        const int m,
                        double *a, const int lda,
                        int *info) {
    F77_FUNC_DPOTRF(&uplo,
                    &m,
                    a, &lda,
                    info);
    return 0;
  }
#if defined(KOKKOS_ENABLE_CUDA)
  template<>
  int 
  Lapack<double>::potrf_buffersize(cusolverDnHandle_t handle,
                                   const cublasFillMode_t uplo,
                                   const int m,
                                   double *a, const int lda,
                                   int *lwork) {
    const int r_val = cusolverDnDpotrf_bufferSize(handle,
                                                  uplo,
                                                  m, 
                                                  a, lda,
                                                  lwork);
    return r_val;
  }

  template<>
  int 
  Lapack<double>::potrf(cusolverDnHandle_t handle,
                        const cublasFillMode_t uplo,
                        const int m,
                        double *a, const int lda,
                        double *w, const int lwork,
                        int *dev) {
    const int r_val = cusolverDnDpotrf(handle,
                                       uplo,
                                       m, 
                                       a, lda,
                                       w, lwork,
                                       dev);
    return r_val;
  }
#endif

  template<>
  int 
  Lapack<double>::sytrf(const char uplo,
                        const int m,
                        double *a, const int lda,
                        int *ipiv,
                        double *work, int lwork,
                        int *info) {
    F77_FUNC_DSYTRF(&uplo,
                    &m,
                    a, &lda,
                    ipiv,
                    work, &lwork,
                    info);
    return 0;
  }

  template<>
  int 
  Lapack<Kokkos::complex<float> >::potrf(const char uplo,
                                         const int m,
                                         Kokkos::complex<float> *a, const int lda,
                                         int *info) {
    F77_FUNC_CPOTRF(&uplo,
                    &m,
                    a, &lda,
                    info);
    return 0;
  }
#if defined(KOKKOS_ENABLE_CUDA)
  template<>
  int 
  Lapack<Kokkos::complex<float> >::potrf_buffersize(cusolverDnHandle_t handle,
                                                    const cublasFillMode_t uplo,
                                                    const int m,
                                                    Kokkos::complex<float> *a, const int lda,
                                                    int *lwork) {
    const int r_val = cusolverDnCpotrf_bufferSize(handle,
                                                  uplo,
                                                  m, 
                                                  (cuComplex*)a, lda,
                                                  lwork);
    return r_val;
  }

  template<>
  int 
  Lapack<Kokkos::complex<float> >::potrf(cusolverDnHandle_t handle,
                                         const cublasFillMode_t uplo,
                                         const int m,
                                         Kokkos::complex<float> *a, const int lda,
                                         Kokkos::complex<float> *w, const int lwork,
                                         int *dev) {
    const int r_val = cusolverDnCpotrf(handle,
                                       uplo,
                                       m, 
                                       (cuComplex*)a, lda,
                                       (cuComplex*)w, lwork,
                                       dev);
    return r_val;
  }
#endif
    
  template<>
  int 
  Lapack<Kokkos::complex<float> >::sytrf(const char uplo,
                                         const int m,
                                         Kokkos::complex<float> *a, const int lda,
                                         int *ipiv,
                                         Kokkos::complex<float> *work, int lwork,
                                         int *info) {
    F77_FUNC_CSYTRF(&uplo,
                    &m,
                    a, &lda,
                    ipiv,
                    work, &lwork,
                    info);
    return 0;
  }

  template<>
  int 
  Lapack<Kokkos::complex<double> >::potrf(const char uplo,
                                          const int m,
                                          Kokkos::complex<double> *a, const int lda,
                                          int *info) {
    F77_FUNC_ZPOTRF(&uplo,
                    &m,
                    a, &lda,
                    info);
    return 0;
  }
#if defined(KOKKOS_ENABLE_CUDA)
  template<>
  int 
  Lapack<Kokkos::complex<double> >::potrf_buffersize(cusolverDnHandle_t handle,
                                                     const cublasFillMode_t uplo,
                                                     const int m,
                                                     Kokkos::complex<double> *a, const int lda,
                                                     int *lwork) {
    const int r_val = cusolverDnZpotrf_bufferSize(handle,
                                                  uplo,
                                                  m, 
                                                  (cuDoubleComplex*)a, lda,
                                                  lwork);
    return r_val;
  }

  template<>
  int 
  Lapack<Kokkos::complex<double> >::potrf(cusolverDnHandle_t handle,
                                          const cublasFillMode_t uplo,
                                          const int m,
                                          Kokkos::complex<double> *a, const int lda,
                                          Kokkos::complex<double> *w, const int lwork,
                                          int *dev) {
    const int r_val = cusolverDnZpotrf(handle,
                                       uplo,
                                       m, 
                                       (cuDoubleComplex*)a, lda,
                                       (cuDoubleComplex*)w, lwork,
                                       dev);
    return r_val;
  }
#endif

  template<>
  int 
  Lapack<Kokkos::complex<double> >::sytrf(const char uplo,
                                          const int m,
                                          Kokkos::complex<double> *a, const int lda,
                                          int *ipiv,
                                          Kokkos::complex<double>* work, int lwork,
                                          int *info) {
    F77_FUNC_ZSYTRF(&uplo,
                    &m,
                    a, &lda,
                    ipiv,
                    work, &lwork,
                    info);
    return 0;
  }

  template<>
  int 
  Lapack<std::complex<float> >::potrf(const char uplo,
                                      const int m,
                                      std::complex<float> *a, const int lda,
                                      int *info) {
    F77_FUNC_CPOTRF(&uplo,
                    &m,
                    (Kokkos::complex<float>*)a, &lda,
                    info);
    return 0;
  }
    
  template<>
  int 
  Lapack<std::complex<float> >::sytrf(const char uplo,
                                      const int m,
                                      std::complex<float> *a, const int lda,
                                      int *ipiv,
                                      std::complex<float> *work, int lwork,
                                      int *info) {
    F77_FUNC_CSYTRF(&uplo,
                    &m,
                    (Kokkos::complex<float>*)a, &lda,
                    ipiv,
                    (Kokkos::complex<float>*)work, &lwork,
                    info);
    return 0;
  }

  template<>
  int 
  Lapack<std::complex<double> >::potrf(const char uplo,
                                       const int m,
                                       std::complex<double> *a, const int lda,
                                       int *info) {
    F77_FUNC_ZPOTRF(&uplo,
                    &m,
                    (Kokkos::complex<double>*)a, &lda,
                    info);
    return 0;
  }
  template<>
  int 
  Lapack<std::complex<double> >::sytrf(const char uplo,
                                       const int m,
                                       std::complex<double> *a, const int lda,
                                       int *ipiv,
                                       std::complex<double>* work, int lwork,
                                       int *info) {
    F77_FUNC_ZSYTRF(&uplo,
                    &m,
                    (Kokkos::complex<double>*)a, &lda,
                    ipiv,
                    (Kokkos::complex<double>*)work, &lwork,
                    info);
    return 0;
  }

}
