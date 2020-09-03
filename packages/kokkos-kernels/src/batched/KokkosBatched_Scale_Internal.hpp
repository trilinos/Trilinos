#ifndef __KOKKOSBATCHED_SCALE_INTERNAL_HPP__
#define __KOKKOSBATCHED_SCALE_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// Serial Internal Impl
  /// ==================== 
  struct SerialScaleInternal {
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const int m, 
           const ScalarType alpha, 
           /* */ ValueType *__restrict__ A, const int as0) {

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i)
        A[i*as0] *= alpha;
        
      return 0;
    }
      
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const int m, const int n, 
           const ScalarType alpha, 
           /* */ ValueType *__restrict__ A, const int as0, const int as1) {

      if (as0 > as1)
        for (int i=0;i<m;++i)
          invoke(n, alpha, A+i*as0, as1);
      else
        for (int j=0;j<n;++j)
          invoke(m, alpha, A+j*as1, as0);
        
      return 0;
    }
  };

  ///
  /// Team Internal Impl
  /// ==================== 
  struct TeamScaleInternal {
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const int m, 
           const ScalarType alpha, 
           /* */ ValueType *__restrict__ A, const int as0) {

      Kokkos::parallel_for
        (Kokkos::TeamThreadRange(member,m),
         [&](const int &i) {
          A[i*as0] *= alpha;
        });
      //member.team_barrier();
      return 0;
    }
      
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const int m, const int n, 
           const ScalarType alpha, 
           /* */ ValueType *__restrict__ A, const int as0, const int as1) {
      if (m > n) {
        Kokkos::parallel_for
          (Kokkos::TeamThreadRange(member,m),
           [&](const int &i) {
            SerialScaleInternal::invoke(n, alpha, A+i*as0, as1);
          });
      } else {
        Kokkos::parallel_for
          (Kokkos::TeamThreadRange(member,n),
           [&](const int &j) {
            SerialScaleInternal::invoke(m, alpha, A+j*as1, as0);
          });
      }
      //member.team_barrier();
      return 0;
    }
  };

  ///
  /// TeamVector Internal Impl
  /// ======================== 
  struct TeamVectorScaleInternal {
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const int m, 
           const ScalarType alpha, 
           /* */ ValueType *__restrict__ A, const int as0) {

      Kokkos::parallel_for
        (Kokkos::TeamVectorRange(member,m),
         [&](const int &i) {
          A[i*as0] *= alpha;
        });
      //member.team_barrier();
      return 0;
    }
      
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const int m, const int n, 
           const ScalarType alpha, 
           /* */ ValueType *__restrict__ A, const int as0, const int as1) {
      if (as0 > as1) {
        Kokkos::parallel_for
          (Kokkos::TeamThreadRange(member,m),
           [&](const int &i) {
            Kokkos::parallel_for
              (Kokkos::ThreadVectorRange(member,n),
               [&](const int &j) {
                A[i*as0+j*as1] *= alpha;
              });
          });
      } else {
        Kokkos::parallel_for
          (Kokkos::ThreadVectorRange(member,m),
           [&](const int &i) {
            Kokkos::parallel_for
              (Kokkos::TeamThreadRange(member,n),
               [&](const int &j) {
                A[i*as0+j*as1] *= alpha;
              });
          });
      }
      //member.team_barrier();
      return 0;
    }
  };

}


#endif
