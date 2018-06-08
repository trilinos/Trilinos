#ifndef __KOKKOSBATCHED_TRSM_TEAM_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRSM_TEAM_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

#include "KokkosBatched_InnerTrsm_Serial_Impl.hpp"
#include "KokkosBatched_Gemm_Team_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {
    ///
    /// Team Internal Impl
    /// ====================

    template<typename AlgoType>
    struct TeamTrsmInternalLeftLower {
      template<typename MemberType,
               typename ScalarType,
               typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const bool use_unit_diag,
             const int m, const int n, 
             const ScalarType alpha,
             const ValueType *__restrict__ A, const int as0, const int as1,
             /**/  ValueType *__restrict__ B, const int bs0, const int bs1);
    };

    template<>
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    TeamTrsmInternalLeftLower<Algo::Trsm::Unblocked>::
    invoke(const MemberType &member, 
           const bool use_unit_diag,
           const int m, const int n,
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ B, const int bs0, const int bs1) {

      const ScalarType one(1.0), zero(0.0);

      if (alpha == zero)   TeamSetInternal  ::invoke(member, m, n, zero,  B, bs0, bs1);
      else {
        if (alpha != one)  TeamScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
        if (m <= 0 || n <= 0) return 0;

        for (int p=0;p<m;++p) {
          const int iend = m-p-1, jend = n;
          
          const ValueType
            *__restrict__ a21 = iend ? A+(p+1)*as0+p*as1 : NULL;
            
          ValueType
            *__restrict__ b1t =        B+p*bs0,
            *__restrict__ B2  = iend ? B+(p+1)*bs0 : NULL;

          member.team_barrier();
          if (!use_unit_diag) {
            const ValueType alpha11 = A[p*as0+p*as1];
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,jend),[&](const int &j) {
                b1t[j*bs1] = b1t[j*bs1] / alpha11;
              });
            member.team_barrier();
          }
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,iend*jend),[&](const int &ij) {
              // assume layout right for batched computation
              const int i = ij/jend, j = ij%jend;
              B2[i*bs0+j*bs1] -= a21[i*as0] * b1t[j*bs1];
            });          
        }
      }      
      return 0;
    }

    template<>
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    TeamTrsmInternalLeftLower<Algo::Trsm::Blocked>::
    invoke(const MemberType &member, 
           const bool use_unit_diag,
           const int m, const int n,
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ B, const int bs0, const int bs1) {

      enum : int {
        mbAlgo = Algo::Trsm::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>()
      };

      const ScalarType one(1.0), zero(0.0), minus_one(-1.0);

      // note that parallel range is different ( m*n vs m-1*n);        
      if (alpha == zero)  TeamSetInternal  ::invoke(member, m, n, zero,  B, bs0, bs1);
      else {
        if (alpha != one) TeamScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
        if (m <= 0 || n <= 0) return 0;

        ///
        /// case host: team size is small and blocksize (mb,nb) is large
            
        ///
        /// case cuda: team size is large and blocksize (mb,nb) is small
        InnerTrsmLeftLowerUnitDiag<mbAlgo>    trsm_u(as0, as1, bs0, bs1);
        InnerTrsmLeftLowerNonUnitDiag<mbAlgo> trsm_n(as0, as1, bs0, bs1);
            
        auto trsm = [&](const int ib, 
                        const int jb,
                        const ValueType *__restrict__ AA,
                        /**/  ValueType *__restrict__ BB) {
          const int mb = mbAlgo;
          const int tsize = member.team_size();
          const int nb = (jb/tsize + jb%tsize > 0);
          const int np = jb%nb;
          for (int p=0;p<ib;p+=mb) {
            const int pb = ((p+mb) > ib ? (ib-p) : mb); 
                
            // trsm update
            const ValueType *__restrict__ Ap = AA+p*as0+p*as1;
            /**/  ValueType *__restrict__ Bp = BB+p*bs0;
                
            member.team_barrier();                  
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,(jb/nb)+(np>0)),[&](const int &jj) {
                const int j = jj*nb, qb = (j+nb) > jb ? np : nb;
                if (use_unit_diag) trsm_u.serial_invoke(Ap, pb, qb, Bp+j*bs1);
                else               trsm_n.serial_invoke(Ap, pb, qb, Bp+j*bs1);
              });
            member.team_barrier();
                
            // gemm update
            TeamGemmInternal<Algo::Gemm::Blocked>
              ::invoke(member,
                       ib-p-pb, jb, pb,
                       minus_one,
                       Ap+pb*as0, as0, as1,
                       Bp, bs0, bs1,
                       one,
                       Bp+pb*bs0, bs0, bs1);
          }
        };
            
        const bool is_small = true; //(m*n <= 64*64);
        if (is_small) {
          trsm(m, n, A, B);
        } else {
          // // some cache blocking may need (not priority yet);
          // trsm(m, n, A, B);
        }
      }
      return 0;
    }

    template<typename AlgoType>
    struct TeamTrsmInternalLeftUpper {
      template<typename MemberType,
               typename ScalarType,
               typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const bool use_unit_diag,
             const int m, const int n, 
             const ScalarType alpha,
             const ValueType *__restrict__ A, const int as0, const int as1,
             /**/  ValueType *__restrict__ B, const int bs0, const int bs1);
    };

    template<>
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::
    invoke(const MemberType &member, 
           const bool use_unit_diag,
           const int m, const int n,
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ B, const int bs0, const int bs1) {

      const ScalarType one(1.0), zero(0.0);

      // note that parallel range is different ( m*n vs m-1*n);        
      if (alpha == zero)  TeamSetInternal  ::invoke(member, m, n, zero,  B, bs0, bs1);
      else {
        if (alpha != one) TeamScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
        if (m <= 0 || n <= 0) return 0;
        
        ValueType *__restrict__ B0 = B;
        for (int p=(m-1);p>=0;--p) {
          const int iend = p, jend = n;

          const ValueType *__restrict__ a01 = A+p*as1;
          /**/  ValueType *__restrict__ b1t = B+p*bs0;
            
          member.team_barrier();
          if (!use_unit_diag) {
            const ValueType alpha11 = A[p*as0+p*as1];
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,jend),[&](const int &j) {
                b1t[j*bs1] = b1t[j*bs1] / alpha11;
              });
            member.team_barrier();
          }

          Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,iend*jend),[&](const int &ij) {
#if							\
  defined (KOKKOS_ENABLE_CUDA) &&				\
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
              const int i = ij%iend, j = ij/iend;
#else
              const int i = ij/jend, j = ij%jend;
#endif
              B0[i*bs0+j*bs1] -= a01[i*as0] * b1t[j*bs1];
            });          
        }
      }
      return 0;
    }

    template<>
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::
    invoke(const MemberType &member,
           const bool use_unit_diag,
           const int m, const int n,
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ B, const int bs0, const int bs1) {

      enum : int {
        mbAlgo = Algo::Trsm::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>()
      };

      const ScalarType one(1.0), zero(0.0), minus_one(-1.0);

      // note that parallel range is different ( m*n vs m-1*n);        
      if (alpha == zero)  TeamSetInternal  ::invoke(member, m, n, zero,  B, bs0, bs1);
      else {
        if (alpha != one) TeamScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
        if (m <= 0 || n <= 0) return 0;

        InnerTrsmLeftUpperUnitDiag<mbAlgo>    trsm_u(as0, as1, bs0, bs1);
        InnerTrsmLeftUpperNonUnitDiag<mbAlgo> trsm_n(as0, as1, bs0, bs1);
            
        auto trsm = [&](const int ib, 
                        const int jb,
                        const ValueType *__restrict__ AA,
                        /**/  ValueType *__restrict__ BB) {
          const int mb = mbAlgo; //(ib <=5 ? ib : mbAlgo);
          const int tsize = member.team_size();
          const int nb = (jb/tsize + jb%tsize > 0);
          const int np = jb%nb;
          for (int pp=0;pp<ib;pp+=mb) {
            const int 
              ptmp = (ib - pp - mb), 
              p = (ptmp < 0 ? 0 : ptmp), 
              pb = (mb + (ptmp < 0)*ptmp);
                  
            // trsm update
            const ValueType *__restrict__ Ap = AA+p*as0+p*as1;
            /**/  ValueType *__restrict__ Bp = BB+p*bs0;

            member.team_barrier();
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,(jb/nb)+(np>0)),[&](const int &jj) {
                const int j = jj*nb, qb = (j+nb) > jb ? np : nb;     
                if (use_unit_diag) trsm_u.serial_invoke(Ap, pb, qb, Bp+j*bs1);
                else               trsm_n.serial_invoke(Ap, pb, qb, Bp+j*bs1);
              });
            member.team_barrier();
                  
            // gemm update
            TeamGemmInternal<Algo::Gemm::Blocked>
              ::invoke(member,
                       p, jb, pb,
                       minus_one,
                       Ap-p*as0, as0, as1,
                       Bp, bs0, bs1,
                       one,
                       BB, bs0, bs1);
          }
        };
          
        const bool is_small = true; //(m*n <= 64*64);
        if (is_small) {
          trsm(m, n, A, B);
        } else {
          // // some cache blocking may need (not priority yet);
          // trsm(m, n, A, B);
        }
      }        
      return 0;      
    }

  }
}

#endif
