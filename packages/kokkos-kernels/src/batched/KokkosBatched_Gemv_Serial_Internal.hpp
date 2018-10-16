#ifndef __KOKKOSBATCHED_GEMV_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMV_SERIAL_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

#include "KokkosBatched_InnerMultipleDotProduct_Serial_Impl.hpp"


namespace KokkosBatched {
  namespace Experimental {

    ///
    /// Serial Internal Impl
    /// ====================

    template<typename ArgAlgo>
    struct SerialGemvInternal {
      template<typename ScalarType,
               typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const int m, const int n, 
             const ScalarType alpha,
             const ValueType *__restrict__ A, const int as0, const int as1,
             const ValueType *__restrict__ x, const int xs0, 
             const ScalarType beta,
             /**/  ValueType *__restrict__ y, const int ys0);
    };

    template<>
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemvInternal<Algo::Gemv::Unblocked>::
    invoke(const int m, const int n, 
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           const ValueType *__restrict__ x, const int xs0,
           const ScalarType beta,
           /**/  ValueType *__restrict__ y, const int ys0) {

      const ScalarType one(1.0), zero(0.0);

      // y = beta y + alpha A x
      // y (m), A(m x n), B(n)

      if      (beta == zero) SerialSetInternal  ::invoke(m, zero, y, ys0);
      else if (beta != one)  SerialScaleInternal::invoke(m, beta, y, ys0);
      
      if (alpha != zero) {
        if (m <= 0 || n <= 0) return 0;

        for (int i=0;i<m;++i) {
          ValueType t(0);
          const ValueType *__restrict__ tA = (A + i*as0);

#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif               
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
          for (int j=0;j<n;++j)
            t += tA[j*as1]*x[j*xs0];
          y[i*ys0] += alpha*t;
        }
      }
      return 0;
    }

    template<>
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemvInternal<Algo::Gemv::Blocked>::
    invoke(const int m, const int n, 
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           const ValueType *__restrict__ x, const int xs0,
           const ScalarType beta,
           /**/  ValueType *__restrict__ y, const int ys0) {

      const ScalarType one(1.0), zero(0.0);

      // y = beta y + alpha A x
      // y (m), A(m x n), B(n)

      enum : int {
        mbAlgo = Algo::Gemv::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>()
      };

      if      (beta == zero) SerialSetInternal  ::invoke(m, zero, y, ys0);
      else if (beta != one)  SerialScaleInternal::invoke(m, beta, y, ys0);
      
      if (alpha != zero) {
        if (m <= 0 || n <= 0) return 0;

        InnerMultipleDotProduct<mbAlgo> inner(as0, as1, xs0, ys0);
        const int mb = mbAlgo;
        for (int i=0;i<m;i+=mb) 
          inner.serial_invoke(alpha, A+i*as0, x, (i+mb) > m ? (m-i) : mb, n, y+i*ys0 );
      }
      return 0;
    }

  }
}


#endif
