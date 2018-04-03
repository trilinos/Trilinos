#ifndef __TACHO_LAPACK_EXTERNAL_HPP__
#define __TACHO_LAPACK_EXTERNAL_HPP__

/// \file  Tacho_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "ShyLU_NodeTacho_config.h"

namespace Tacho {

    template<typename T>
    struct Lapack {
      static 
      void potrf(const char uplo,
                 const int m,
                 T *a, const int lda,
                 int *info);

      static 
      void sytrf(const char uplo,
                 const int m,
                 T *a, const int lda,
                 int *ipiv,
                 T *work, int lwork,
                 int *info);
    };
}

#endif
