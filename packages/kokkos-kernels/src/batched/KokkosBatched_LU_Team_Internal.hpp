#ifndef __KOKKOSBATCHED_LU_TEAM_INTERNAL_HPP__
#define __KOKKOSBATCHED_LU_TEAM_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_InnerLU_Serial_Impl.hpp"
#include "KokkosBatched_InnerTrsm_Serial_Impl.hpp"

#include "KokkosBatched_Trsm_Team_Internal.hpp"
#include "KokkosBatched_Gemm_Team_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {
  
    ///
    /// Team Internal Impl
    /// ==================

    template<typename AlgoType>
    struct TeamLU_Internal {
      template<typename MemberType, typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static int 
      invoke(const MemberType &member,
             const int m, const int n,
             ValueType *__restrict__ A, const int as0, const int as1);
    };

    template<>
    template<typename MemberType, typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    TeamLU_Internal<Algo::LU::Unblocked>::
    invoke(const MemberType &member, 
           const int m, const int n,
           ValueType *__restrict__ A, const int as0, const int as1) {

      const int k = (m < n ? m : n);
      if (k <= 0) return 0;

      for (int p=0;p<k;++p) {
        const int iend = m-p-1, jend = n-p-1;

        const ValueType 
          // inv_alpha11 = 1.0/A(p,p),
          alpha11 = A[p*as0+p*as1],
          *__restrict__ a12t = A+(p  )*as0+(p+1)*as1;
            
        ValueType
          *__restrict__ a21  = A+(p+1)*as0+(p  )*as1,
          *__restrict__ A22  = A+(p+1)*as0+(p+1)*as1;
            
        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,iend),[&](const int &i) {
            // a21[i*as0] *= inv_alpha11; 
            a21[i*as0] /= alpha11;
          });
            
        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,iend*jend),[&](const int &ij) {
#if							\
  defined (KOKKOS_HAVE_CUDA) &&				\
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
            const int i = ij%iend, j = ij/iend;
#else
            const int i = ij/jend, j = ij%jend;
#endif
            A22[i*as0+j*as1] -= a21[i*as0] * a12t[j*as1];
          });
      }
      return 0;
    }
    
    template<>
    template<typename MemberType, typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    TeamLU_Internal<Algo::LU::Blocked>::
    invoke(const MemberType &member, 
           const int m, const int n,
           ValueType *__restrict__ A, const int as0, const int as1) {

      enum : int {
        mbAlgo = Algo::LU::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>()
      };

      const int k = (m < n ? m : n);
      if (k <= 0) return 0;

      const typename Kokkos::Details::ArithTraits<ValueType>::mag_type one(1.0), minus_one(-1.0);

      InnerLU<mbAlgo> lu(as0, as1);
          
      InnerTrsmLeftLowerUnitDiag<mbAlgo>    trsm_llu(as0, as1, as0, as1);
      InnerTrsmLeftLowerNonUnitDiag<mbAlgo> trsm_run(as1, as0, as1, as0);
      auto lu_factorize = [&](const int ib,
                              const int jb,
                              ValueType *__restrict__ AA) {
        const int tsize = member.team_size();
        const int mb = mbAlgo;
        const int nb = ((jb-mb) + (ib-mb))/tsize + ((jb-mb) + (ib-mb))%tsize > 0;
        const int kb = ib < jb ? ib : jb; 

        for (int p=0;p<kb;p+=mb) {
          const int pb = (p+mb) > kb ? (kb-p) : mb;

          // diagonal block
          ValueType *__restrict__ Ap = AA+p*as0+p*as1;
                
          // lu on a block             
          member.team_barrier();
          if (member.team_rank() == 0)
            lu.serial_invoke(pb, Ap);
          member.team_barrier();

          const int 
            m_abr  = ib-p-mb,               n_abr  = jb-p-mb,
            mp_abr = m_abr%nb,              np_abr = n_abr%nb,
            mq_abr = (m_abr/nb)+(mp_abr>0), nq_abr = (n_abr/nb)+(np_abr>0);

          // trsm update
          Kokkos::parallel_for
            (Kokkos::TeamThreadRange(member,0,mq_abr+nq_abr),
             [&](const int &ij) {
              if (ij < nq_abr) {
                const int j = (ij)*nb, qb = (j+nb) > n_abr ? np_abr : nb;
                trsm_llu.serial_invoke(Ap, pb, qb, Ap+(j+mb)*as1);
              } else {
                const int i = (ij-nq_abr)*nb , qb = (i+nb) > m_abr ? mp_abr : nb;
                trsm_run.serial_invoke(Ap, pb, qb, Ap+(i+mb)*as0);
              }
            });
          member.team_barrier();

          // gemm update
          TeamGemmInternal<Algo::Gemm::Blocked>::
            invoke(member, 
                   m_abr, n_abr, pb,
                   minus_one,
                   Ap+mb*as0, as0, as1,
                   Ap+mb*as1, as0, as1,
                   one,
                   Ap+mb*as0+mb*as1, as0, as1);
        }
      };
            
      const bool is_small = true; //(m*n <= 64*64);
      if (is_small) {
        lu_factorize(m, n, A);
      } else {
        // // some cache blocking may need (not priority yet);
        // lu_factorize(m, n, A);
      }

      return 0;
    }
  }
}


#endif
