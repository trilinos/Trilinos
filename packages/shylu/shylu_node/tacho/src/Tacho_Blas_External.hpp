#ifndef __TACHO_BLAS_EXTERNAL_HPP__
#define __TACHO_BLAS_EXTERNAL_HPP__

/// \file  Tacho_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "ShyLU_NodeTacho_config.h"

namespace Tacho {

    template<typename T>
    struct Blas {
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
      void trsm(const char side, const char uplo, const char transa, const char diag,
                int m, int n, 
                const T alpha, 
                const T *a, int lda,
                /* */ T *b, int ldb);
    };

}

#endif
