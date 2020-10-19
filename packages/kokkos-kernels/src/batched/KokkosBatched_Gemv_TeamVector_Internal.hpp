#ifndef __KOKKOSBATCHED_GEMV_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMV_TEAMVECTOR_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

#include "KokkosBatched_InnerMultipleDotProduct_Serial_Impl.hpp"

namespace KokkosBatched {

  ///
  /// Team Internal Impl
  /// ====================
  template<typename ArgAlgo>
  struct TeamVectorGemvInternal {
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const int m, const int n, 
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           const ValueType *__restrict__ x, const int xs0, 
           const ScalarType beta,
           /**/  ValueType *__restrict__ y, const int ys0) {
      assert(false && "Error: encounter dummy impl");
      return 0;
    }
  };
    
  template<>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  TeamVectorGemvInternal<Algo::Gemv::Unblocked>::
  invoke(const MemberType &member,
         const int m, const int n, 
         const ScalarType alpha,
         const ValueType *__restrict__ A, const int as0, const int as1,
         const ValueType *__restrict__ x, const int xs0,
         const ScalarType beta,
         /**/  ValueType *__restrict__ y, const int ys0) {
    const ScalarType one(1.0), zero(0.0);

    // y = beta y + alpha A x
    // y (m), A(m x n), B(n)
      
    if      (beta == zero) TeamVectorSetInternal  ::invoke(member, m, zero, y, ys0);
    else if (beta != one)  TeamVectorScaleInternal::invoke(member, m, beta, y, ys0);
      
    if (alpha != zero) {
      if (m <= 0 || n <= 0) return 0;
        
      if (beta != one) 
        member.team_barrier();
        
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,m),[&](const int &i) {
          ValueType t(0);
          const ValueType *__restrict__ tA = (A + i*as0);
          Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(member,n),[&](const int &j, ValueType &update) {
              update += tA[j*as1]*x[j*xs0];
            }, t);
          Kokkos::single(Kokkos::PerThread(member), [&]() {
              y[i*ys0] += alpha*t;
            });
        });
    }
    return 0;
  }

}


#endif
