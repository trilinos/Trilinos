/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels: Linear Algebra and Graph Kernels
//                 Copyright 2016 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/





#ifndef __KOKKOSKERNELS_GEMM_TEAM_IMPL_HPP__
#define __KOKKOSKERNELS_GEMM_TEAM_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_Util.hpp"

namespace KokkosKernels {

  namespace Team {
    template<typename MemberType>
    struct Gemm<MemberType,Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Triple> {

      template<typename ScalarType,
               typename AViewType,
               typename BViewType,
               typename CViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member,
             const ScalarType alpha,
             const AViewType A,
             const BViewType B,
             const ScalarType beta,
             /**/  CViewType C) {
        // C = beta C + alpha A B
        // C (m x n), A(m x k), B(k x n)
        
        typedef typename CViewType::value_type value_type;
        
        const int 
          team_rank = member.team_rank();
        
        // serial update on C
        if (team_rank == 0) {
          if      (beta == 0) Util::set  (C, value_type(0)   );
          else if (beta != 1) Util::scale(C, value_type(beta));
        }
        
        member.team_barrier();
        
        if (alpha != 0) {
          const int
            m = C.dimension(0),
            n = C.dimension(1),
            k = A.dimension(1);
          
          const int
            as0 = A.stride_0(),
            bs1 = B.stride_1(),
            cs0 = C.stride_0(),
            cs1 = C.stride_1();
          
          value_type
            *__restrict__ pC = &C(0,0);
          for (int p=0;p<k;++p) {
            const value_type
              *__restrict__ pA = &A(0,p),
              *__restrict__ pB = &B(p,0);
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,m),[&](const int &i) {
                // for (int i=0;i<m;++i) {
                const value_type tA(alpha*pA[i*as0]);
                //Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,n),[&](const int &j) {
                for (int j=0;j<n;++j)
                  pC[i*cs0+j*cs1] += tA*pB[j*bs1];
              });
          }
        }
        return 0;
      }
    };
    
    template<typename MemberType>
    struct Gemm<MemberType,Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType,
               typename CViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const ScalarType alpha,
             const AViewType A,
             const BViewType B,
             const ScalarType beta,
             /**/  CViewType C) {
        // C = beta C + alpha A B
        // C (m x n), A(m x k), B(k x n)
        
        typedef typename CViewType::value_type value_type;
        
        const int 
          team_rank = member.team_rank();
        
        if (team_rank == 0) {
          if      (beta == 0) Util::set  (C, value_type(0)   );
          else if (beta != 1) Util::scale(C, value_type(beta));
        }
        
        if (alpha != 0) {
          const int
            m = C.dimension(0),
            n = C.dimension(1),
            k = A.dimension(1);
          
          const int
            as0 = A.stride_0(),
            as1 = A.stride_1(),
            bs0 = B.stride_0(),
            bs1 = B.stride_1(),
            cs0 = C.stride_0(),
            cs1 = C.stride_1();
          
	  enum : int {
            mb = Algo::Gemm::Blocked::mb,
            nb = Algo::Gemm::Blocked::nb };
          
          InnerRankUpdate<mb,nb> inner(as0, as1,
                                       bs0, bs1,
                                       cs0, bs1);
          const int mm = (m/mb)*mb, nn = (n/nb)*nb;
	  enum : int {
            team_threshold_size = 32 };

          // // no benefit on knl
          // // 0th loop parallel team
          // for (int i=0;i<mm;i+=mb)
          //   for (int j=0;j<nn;j+=nb)
          //     inner.team_invoke_var2(member, alpha, &A(i,0), &B(0,j), k, &C(i,j));

          // small matrix; there is no team efficiency
          if (m <= team_threshold_size && 
              n <= team_threshold_size) {
            // no team parallel
            if (team_rank == 0) 
              for (int i=0;i<mm;i+=mb)
                for (int j=0;j<nn;j+=nb)
                  inner.serial_invoke(alpha, &A(i,0), &B(0,j), k, &C(i,j));
          } 
          else if (m < n) {
            // 1st loop parallel team
            for (int i=0;i<mm;i+=mb)
              Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,nn/nb),[&](const int jj) {
                  const int j = jj*nb;
                  inner.serial_invoke(alpha, &A(i,0), &B(0,j), k, &C(i,j));
                });
          } 
          else {
            // 2nd loop parallel team
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,mm/mb),[&](const int ii) {
                const int i = ii*mb;
                for (int j=0;j<nn;j+=nb)
                  inner.serial_invoke(alpha, &A(i,0), &B(0,j), k, &C(i,j));
              });
          }
          
          if (team_rank == 0) {
            const int mp = (m%mb), np = (n%nb);
            if (mp      ) inner.serial_invoke(alpha, &A(mm, 0), &B( 0, 0), mp, nn, k, &C(mm, 0));
            if (      np) inner.serial_invoke(alpha, &A( 0, 0), &B( 0,nn), mm, np, k, &C( 0,nn));
            if (mp && np) inner.serial_invoke(alpha, &A(mm, 0), &B( 0,nn), mp, np, k, &C(mm,nn));
          }

          member.team_barrier();
        }
        return 0;
      }
    };

  } // end namespace Team

  template<int mb, int nb>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerRankUpdate<mb,nb>::
  team_invoke(const MemberType &member,
              const ScalarType alpha,
              const ValueType *__restrict__ A,
              const ValueType *__restrict__ B,
              const int m, const int n, const int k,
              /**/  ValueType *__restrict__ C) {
    for (int p=0;p<k;++p) {
      const ValueType
        *__restrict__ pA = A + p*_as1,
        *__restrict__ pB = B + p*_bs0;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,m),[&](const int &i) {
          //for (int i=0;i<m;++i) {
          const ValueType tA(alpha*pA[i*_as0]);
          //Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,n),[&](const int &j) {
          for (int j=0;j<n;++j)
            C[i*_cs0+j*_cs1] += tA*pB[j*_bs1];
        });
    }
    return 0;
  }
  
  template<>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerRankUpdate<4,4>::
  team_invoke_var1(const MemberType &member, 
                   const ScalarType alpha,
                   const ValueType *__restrict__ A,
                   const ValueType *__restrict__ B,
                   const int k,
                   /**/  ValueType *__restrict__ C) {
    const int
      i0 = 0*_as0, i1 = 1*_as0, i2 = 2*_as0, i3 = 3*_as0,
      j0 = 0*_bs1, j1 = 1*_bs1, j2 = 2*_bs1, j3 = 3*_bs1;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,4),[&](const int &j) {
        ValueType 
          c_0 = 0,
          c_1 = 0,
          c_2 = 0,
          c_3 = 0;
        
        for (int p=0;p<k;++p) {
          const ValueType 
            a_0 = A[i0+p*_as1], b_0 = B[p*_bs0+j0],
            a_1 = A[i1+p*_as1], b_1 = B[p*_bs0+j1],
            a_2 = A[i2+p*_as1], b_2 = B[p*_bs0+j2],
            a_3 = A[i3+p*_as1], b_3 = B[p*_bs0+j3];
          
          c_0 += a_0 * b_0;
          c_1 += a_1 * b_1;
          c_2 += a_2 * b_2;
          c_3 += a_3 * b_3;
        }
        C[0*_cs0+j*_cs1] += alpha * c_0; 
        C[1*_cs0+j*_cs1] += alpha * c_1; 
        C[2*_cs0+j*_cs1] += alpha * c_2; 
        C[3*_cs0+j*_cs1] += alpha * c_3;      
      });
    
    return 0;
  }

  template<>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerRankUpdate<4,4>::
  team_invoke_var2(const MemberType &member, 
                   const ScalarType alpha,
                   const ValueType *__restrict__ A,
                   const ValueType *__restrict__ B,
                   const int k,
                   /**/  ValueType *__restrict__ C) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,16),[&](const int &ij) {
        const int 
          i = ij/4,
          j = ij%4;
        
        ValueType 
          c = 0;
        
        for (int p=0;p<k;++p) {
          const ValueType 
            a = A[i+p*_as1], b = B[p*_bs0+j];
          c += a * b;
        }
        C[i*_cs0+j*_cs1] += alpha * c; 
      });
    
    return 0;
  }


}

#endif
