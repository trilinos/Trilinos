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
    static 
    int sytrf(const char uplo,
              const int m,
              T *a, const int lda,
              int *ipiv,
              T *work, int lwork,
              int *info);
  };
}

#endif
