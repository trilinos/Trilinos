#ifndef __TACHO_LAPACK_EXTERNAL_HPP__
#define __TACHO_LAPACK_EXTERNAL_HPP__

/// \file  Tacho_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp" // CUDA specialization

#if defined(KOKKOS_ENABLE_CUDA) 
#include "cusolverDn.h"
#endif

namespace Tacho {

  template<typename T>
  struct Lapack {
    ///
    /// Cholesky
    ///
    static 
    int potrf(const char uplo,
              const int m,
              T *a, const int lda,
              int *info);
#if defined (KOKKOS_ENABLE_CUDA)
    static
    int potrf_buffersize(cusolverDnHandle_t handle,
                         const cublasFillMode_t uplo,
                         const int m,
                         T *a, const int lda,
                         int *lwork);
    static
    int potrf(cusolverDnHandle_t handle,
              const cublasFillMode_t uplo,
              const int m,
              T *a, const int lda,
              T *W, const int lwork,
              int *dev);
#endif

    ///
    /// LDLt 
    ///
    static 
    int sytrf(const char uplo,
              const int m,
              T *a, const int lda,
              int *ipiv,
              T *work, int lwork,
              int *info);
#if defined (KOKKOS_ENABLE_CUDA)
    static
    int sytrf_buffersize(cusolverDnHandle_t handle,
                         const int m,
                         T *a, const int lda,
                         int *lwork);
    static
    int sytrf(cusolverDnHandle_t handle,
              const cublasFillMode_t uplo,
              const int m,
              T *a, const int lda,
              int *ipiv,
              T *W, const int lwork,
              int *dev);
#endif

    ///
    /// LU
    ///
    static
    int getrf(const int m, const int n,
              T *a, const int lda,
              int *ipiv,
              int *info);
#if defined (KOKKOS_ENABLE_CUDA)
    static
    int getrf_buffersize(cusolverDnHandle_t handle,
                         const int m, const int n,
                         T *a, const int lda,
                         int *lwork);
    static
    int getrf(cusolverDnHandle_t handle,
              const int m, const int n,
              T *a, const int lda,
              T *w,
              int *ipiv,
              int *dev);
#endif
    
  };
}

#endif
