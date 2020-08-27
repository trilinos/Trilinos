#ifndef __KOKKOSBATCHED_DOT_INTERNAL_HPP__
#define __KOKKOSBATCHED_DOT_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// Serial Internal Impl
  /// ==================== 
        
  struct SerialDotInternal {

    // i \in [0,m)  
    // C = conj(A(:))*B(:)  
    template<typename ValueType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const int m, 
           const ValueType *__restrict__ A, const int as0,
	   const ValueType *__restrict__ B, const int bs0, 
           /* */ ValueType *__restrict__ C) {
      using ats = Kokkos::ArithTraits<ValueType>;
      C[0] = ValueType(0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i) {
	const int idx_a = i*as0, idx_b = i*bs0;
	C[0] += ats::conj(A[idx_a])*B[idx_b];
      }
      return 0;
    }

    // j \in [0,n), i \in [0,m)
    // C(j) = conj(A(:,j))*B(:,j)
    template<typename ValueType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const int m, const int n, 
           const ValueType *__restrict__ A, const int as0, const int as1,
	   const ValueType *__restrict__ B, const int bs0, const int bs1,
           /* */ ValueType *__restrict__ C, const int cs) {
      for (int j=0;j<n;++j)               
	invoke(m, A+j*as1, as0, B+j*bs1, bs0, C+j*cs);
      return 0;
    }
  };        
    
  ///
  /// TeamVector Internal Impl
  /// ========================

  // i \in [0,m)  
  // C = conj(A(:))*B(:)  
  struct TeamVectorDotInternal {
    template<typename MemberType,
             typename ValueType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const int m, 
           const ValueType *__restrict__ A, const int as0,
	   const ValueType *__restrict__ B, const int bs0, 
           /* */ ValueType *__restrict__ C) {
      using ats = Kokkos::ArithTraits<ValueType>;
      ValueType t(0);
      Kokkos::parallel_reduce
        (Kokkos::TeamVectorRange(member,m),
	 [&](const int &i, ValueType &update) {
	   const int idx_a = i*as0, idx_b = i*bs0; 
	   update += ats::conj(A[idx_a])*B[idx_b];
	 }, t);
       Kokkos::single
	 (Kokkos::PerThread(member),
	  [&]() {
	    C[0] = t;
	  });
      return 0;
    }

    // j \in [0,n), i \in [0,m)
    // C(j) = conj(A(:,j))*B(:,j)
    template<typename MemberType,
             typename ValueType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const int m, const int n, 
           const ValueType *__restrict__ A, const int as0, const int as1,
	   const ValueType *__restrict__ B, const int bs0, const int bs1,
           /* */ ValueType *__restrict__ C, const int cs) {
      using ats = Kokkos::ArithTraits<ValueType>;
      Kokkos::parallel_for
	(Kokkos::TeamThreadRange(member,n),
	 [&](const int &j) {
	   ValueType t(0);
	   const ValueType *__restrict__ A_at_j = A + j*as1;
	   const ValueType *__restrict__ B_at_j = B + j*bs1;
	   Kokkos::parallel_reduce
	     (Kokkos::ThreadVectorRange(member,m),
	      [&](const int &i, ValueType &update) {
		const int idx_a = i*as0, idx_b = i*bs0;
		update += ats::conj(A_at_j[idx_a])*B_at_j[idx_b];
	      }, t);
	   Kokkos::single
	     (Kokkos::PerThread(member),
	      [&]() {
		C[j*cs] = t;
	      });
	 });
      return 0;
    }
  };

} // end namespace KokkosBatched


#endif
