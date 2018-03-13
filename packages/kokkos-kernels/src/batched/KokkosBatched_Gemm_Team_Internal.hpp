#ifndef __KOKKOSBATCHED_GEMM_TEAM_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMM_TEAM_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

#include "KokkosBatched_InnerGemmFixC_Serial_Impl.hpp"


namespace KokkosBatched {
  namespace Experimental {
    ///
    /// Team Internal Impl
    /// ==================== 
    template<typename ArgAlgo>
    struct TeamGemmInternal {
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
    TeamGemmInternal<Algo::Gemm::Unblocked>::
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
        
      if      (beta == zero) TeamSetInternal  ::invoke(member, m, n, zero, C, cs0, cs1);
      else if (beta != one ) TeamScaleInternal::invoke(member, m, n, beta, C, cs0, cs1);
        
      if (alpha != ScalarType(0.0)) {
        if (m <= 0 || n <= 0 || k <= 0) return 0;

        if (beta != one) 
          member.team_barrier();
            
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,m*n),[&](const int &ij) {
#if                                                     \
  defined (KOKKOS_ENABLE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
            const int i = ij%m, j = ij/m;
#else
            const int i = ij/n, j = ij%n;
#endif
            const ValueType
              *__restrict__ pA = A+i*as0,
              *__restrict__ pB = B+j*bs1;
            
            ValueType c = 0;
            for (int p=0;p<k;++p) 
              c += pA[p*as1]*pB[p*bs0];
            C[i*cs0+j*cs1] += alpha*c;
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
    TeamGemmInternal<Algo::Gemm::Blocked>::
    invoke(const MemberType &member, 
           const int m, const int n, const int k,
           const ScalarType alpha, 
           const ValueType *__restrict__ A, const int as0, const int as1,
           const ValueType *__restrict__ B, const int bs0, const int bs1,
           const ScalarType beta,
           /**/  ValueType *__restrict__ C, const int cs0, const int cs1) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)

      enum : int {
        mbAlgo = Algo::Gemm::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>(),
        nbAlgo = Algo::Gemm::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>() 
      };

      const ScalarType one(1.0), zero(0.0);
          
      if      (beta == zero) TeamSetInternal  ::invoke(member, m, n, zero, C, cs0, cs1);
      else if (beta != one ) TeamScaleInternal::invoke(member, m, n, beta, C, cs0, cs1);

      if (alpha != ScalarType(0.0)) {
        if (m <= 0 || n <= 0 || k <= 0) return 0;

        if (beta != one)
          member.team_barrier();

        ///
        /// case cuda: team size is large and blocksize (mb,nb) is small
        InnerGemmFixC<mbAlgo,nbAlgo> inner(as0, as1, bs0, bs1, cs0, cs1);
        auto gemm = [&](const int ib, 
                        const int jb,
                        const int pb,
                        const ValueType *__restrict__ AA,
                        const ValueType *__restrict__ BB,
                        /**/  ValueType *__restrict__ CC) {
          const int              
          mb = mbAlgo, mp = (ib%mb), mq = (ib/mb) + (mp>0),
          nb = nbAlgo, np = (jb%nb), nq = (jb/nb) + (np>0);
              
          // square tiling
          Kokkos::parallel_for
          (Kokkos::TeamThreadRange(member, mq*nq ),
           [&](const int &ij) {
#if                                                     \
  defined (KOKKOS_ENABLE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
            const int i = ij%mq*mb, j = ij/mq*nb;
#else
            const int i = ij/nq*mb, j = ij%nq*nb;
#endif
            inner.serial_invoke(alpha, 
                                AA+i*as0, BB+j*bs1, 
                                (i+mb) > ib ? mp : mb, 
                                (j+nb) > jb ? np : nb, 
                                pb, 
                                CC+i*cs0+j*cs1);
          });
        };          
            
        const bool is_small = true; //(m*n*k <= 64*64*64);
        if (is_small) {
          gemm(m, n, k, A, B, C);
        } else {
          // // cache blocking
          // const int 
          //   nc = nb*10, kc = mb*4, mc = mb*4;
          
          // for (int jj=0;jj<n;jj+=nc) {
          //   const int tj = n-jj, jb = (tj < nc ? tj : nc);
          //   for (int pp=0;pp<k;pp+=kc) {
          //     const int tp = k-pp, pb = (tp < kc ? tp : kc);
          //     //const int pb = k, pp = 0;
          //     for (int ii=0;ii<m;ii+=mc) {
          //       const int ti = m-ii, ib = (ti < mc ? ti : mc);
                
          //       const ValueType *__restrict__ AA = A+ii*as0+pp*as1;
          //       const ValueType *__restrict__ BB = B+pp*bs0+jj*bs1;
          //       /**/  ValueType *__restrict__ CC = C+ii*cs0+jj*cs1;
                
          //       gemm(ib, jb, pb, AA, BB, CC);                  
          //     } // for ii
          //   } // for pp
          // } // for jj          
        }

      }
      return 0;
    }

  }
}

#endif
