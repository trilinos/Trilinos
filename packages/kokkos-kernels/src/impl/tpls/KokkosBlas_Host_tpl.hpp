#ifndef KOKKOSBLAS_HOST_TPL_HPP_
#define KOKKOSBLAS_HOST_TPL_HPP_

/// \file  KokkosBlas_Host_tpl.hpp 
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)


#include "KokkosKernels_config.h"
#include "Kokkos_ArithTraits.hpp"

#if defined( KOKKOSKERNELS_ENABLE_TPL_BLAS )

namespace KokkosBlas {
  namespace Impl {

    template<typename T>
    struct HostBlas {
      typedef Kokkos::ArithTraits<T> ats;
      typedef typename ats::mag_type mag_type;
  
      static
      void scal(int n,
                const T alpha,
                /* */ T *x, int x_inc);

      static 
      int iamax(int n, 
                const T *x, int x_inc);
      
      static 
      mag_type nrm2(int n,
                    const T *x, int x_inc);

      static
      mag_type asum(int n, 
                    const T *x, int x_inc);
      
      static 
      T dot(int n,
            const T *x, int x_inc,
            const T *y, int y_inc);

      static
      void axpy(int n, 
                const T alpha,
                const T *x, int x_inc,
                /* */ T *y, int y_inc);

      static 
      void gemv(const char trans, 
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc);

      static 
      void trsv(const char uplo, const char transa, const char diag, 
                int m, 
                const T *a, int lda,
                /* */ T *b, int ldb);

      static 
      void gemm(const char transa, const char transb, 
                int m, int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T *b, int ldb,
                const T beta,
                /* */ T *c, int ldc);

      static 
      void herk(const char transa, const char transb, 
                int n, int k,
                const T alpha, 
                const T *a, int lda,
                const T beta,
                /* */ T *c, int ldc);

      static 
      void trmm(const char side, const char uplo, const char transa, const char diag,
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                /* */ T *b, int ldb);

      static 
      void trsm(const char side, const char uplo, const char transa, const char diag,
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                /* */ T *b, int ldb);

      static
      void gesv(int n, int rhs,
                T *a, int lda, int *ipiv,
                T *b, int ldb,
                int info);
    };

  }
}

#endif

#endif
