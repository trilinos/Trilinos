#ifndef __KOKKOSKERNELS_GEMM_SERIAL_DECL_HPP__
#define __KOKKOSKERNELS_GEMM_SERIAL_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <immintrin.h>

namespace KokkosKernels {

  namespace Serial {

    template<typename ArgTransA,
             typename ArgTransB,
             typename ArgAlgo>
    struct Gemm {

      template<typename ScalarType,
               typename AViewType,
               typename BViewType,
               typename CViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType A,
             const BViewType B,
             const ScalarType beta,
             /**/  CViewType C,
             /**/  void *aux = NULL);
      
    };

  }

}

#endif
