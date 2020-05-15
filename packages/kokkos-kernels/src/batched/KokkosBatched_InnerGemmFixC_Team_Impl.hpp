#ifndef __KOKKOSBATCHED_INNER_GEMM_FIX_C_TEAM_IMPL_HPP__
#define __KOKKOSBATCHED_INNER_GEMM_FIX_C_TEAM_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_InnerGemmFixC_Decl.hpp"

namespace KokkosBatched {

  template<int mb, int nb>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerGemmFixC<mb,nb>::
  team_invoke(const MemberType &member,
              const ScalarType alpha,
              const ValueType *__restrict__ A,
              const ValueType *__restrict__ B,
              const int k,
              /**/  ValueType *__restrict__ C) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,mb*nb),[&](const int &ij) {
        const int
          i = ij/nb,
          j = ij%nb;

        const ValueType
          *__restrict__ pA = A+i*_as0,
          *__restrict__ pB = B+j*_bs1;

        ValueType c = 0;
        for (int p=0;p<k;++p)
          c += pA[p*_as1]*pB[p*_bs0];
        C[i*_cs0+j*_cs1] += alpha*c;
      });
    return 0;
  }

  template<int mb, int nb>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerGemmFixC<mb,nb>::
  team_invoke(const MemberType &member,
              const ScalarType alpha,
              const ValueType *__restrict__ A,
              const ValueType *__restrict__ B,
              const int m, const int n, const int k,
              /**/  ValueType *__restrict__ C) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,m*n),[&](const int &ij) {
        const int
          i = ij/n,
          j = ij%n;

        const ValueType
          *__restrict__ pA = A+i*_as0,
          *__restrict__ pB = B+j*_bs1;

        ValueType c = 0;
        for (int p=0;p<k;++p)
          c += pA[p*_as1]*pB[p*_bs0];
        C[i*_cs0+j*_cs1] += alpha*c;
      });
    return 0;
  }
}


#endif
