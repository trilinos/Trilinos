#ifndef __KOKKOSBATCHED_COPY_INTERNAL_HPP__
#define __KOKKOSBATCHED_COPY_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// Serial Internal Impl
  /// ==================== 
        
  struct SerialCopyInternal {
    template<typename ValueType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const int m, 
           const ValueType *__restrict__ A, const int as0, 
           /* */ ValueType *__restrict__ B, const int bs0) {
        
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i) 
        B[i*bs0] = A[i*as0];
        
      return 0;
    }
    template<typename ValueType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const int m, const int n, 
           const ValueType *__restrict__ A, const int as0, const int as1,
           /* */ ValueType *__restrict__ B, const int bs0, const int bs1) {
      if (as1 < as0) 
        for (int i=0;i<m;++i) 
          invoke(n, A+i*as0, as1, B+i*bs0, bs1);
      else 
        for (int j=0;j<n;++j)               
          invoke(m, A+j*as1, as0, B+j*bs1, bs0);
      return 0;
    }
  };        
    
  ///
  /// Team Internal Impl
  /// ==================
  struct TeamCopyInternal {
    template<typename MemberType,
             typename ValueType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const int m, 
           const ValueType *__restrict__ A, const int as0, 
           /* */ ValueType *__restrict__ B, const int bs0) {
      Kokkos::parallel_for
        (Kokkos::TeamThreadRange(member,0,m),[&](const int &i) {
          B[i*bs0] = A[i*as0];
        });
      //member.team_barrier();
      return 0;
    }
    template<typename MemberType,
             typename ValueType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const int m, const int n, 
           const ValueType *__restrict__ A, const int as0, const int as1,
           /* */ ValueType *__restrict__ B, const int bs0, const int bs1) {
      if (m > n) { 
        Kokkos::parallel_for
          (Kokkos::TeamThreadRange(member,0,m),[&](const int &i) {
            SerialCopyInternal::invoke(n, A+i*as0, as1, B+i*bs0, bs1);
          });
      } else {
        Kokkos::parallel_for
          (Kokkos::TeamThreadRange(member,0,n),[&](const int &j) {
            SerialCopyInternal::invoke(m, A+j*as1, as0, B+j*bs1, bs0);
          });
      }
      //member.team_barrier();
      return 0;
    }
  };

} // end namespace KokkosBatched


#endif
