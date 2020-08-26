#ifndef __KOKKOSBATCHED_GEMM_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMM_TEAMVECTOR_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

namespace KokkosBatched {

  ///
  /// TeamVector Internal Impl
  /// ==================== 
  template<typename ArgAlgo>
  struct TeamVectorGemmInternal {
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const int m, const int n, const int k,
           const ScalarType alpha, 
           const ValueType *__restrict__ A, const int as0, const int as1,
           const ValueType *__restrict__ B, const int bs0, const int bs1,
           const ScalarType beta,
           /**/  ValueType *__restrict__ C, const int cs0, const int cs1);
  };

  template<>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  TeamVectorGemmInternal<Algo::Gemm::Unblocked>::
  invoke(const MemberType &member, 
         const int m, const int n, const int k,
         const ScalarType alpha, 
         const ValueType *__restrict__ A, const int as0, const int as1,
         const ValueType *__restrict__ B, const int bs0, const int bs1,
         const ScalarType beta,
         /**/  ValueType *__restrict__ C, const int cs0, const int cs1) {

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
      
    const ScalarType one(1.0), zero(0.0);
        
    if      (beta == zero) TeamVectorSetInternal  ::invoke(member, m, n, zero, C, cs0, cs1);
    else if (beta != one ) TeamVectorScaleInternal::invoke(member, m, n, beta, C, cs0, cs1);
        
    if (alpha != ScalarType(0.0)) {
      if (m <= 0 || n <= 0 || k <= 0) return 0;

      if (beta != one) 
        member.team_barrier();
            
      Kokkos::parallel_for
	(Kokkos::TeamThreadRange(member,m),
	 [&](const int &i) {
	   const ValueType
	     *__restrict__ pA = A+i*as0;
	   Kokkos::parallel_for
	     (Kokkos::ThreadVectorRange(member,n),
	      [&](const int &j) {								   
		const ValueType
		  *__restrict__ pB = B+j*bs1;
		
		ValueType c = 0;
		for (int p=0;p<k;++p) 
		  c += pA[p*as1]*pB[p*bs0];
		C[i*cs0+j*cs1] += alpha*c;
	      });
	 });
    }
    return 0;
  }

}


#endif
