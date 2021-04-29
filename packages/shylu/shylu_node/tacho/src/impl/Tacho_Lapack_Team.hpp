#ifndef __TACHO_LAPACK_TEAM_HPP__
#define __TACHO_LAPACK_TEAM_HPP__

/// \file  Tacho_Lapack_TEAM.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp"

namespace Tacho {

    template<typename T>
    struct LapackTeam {
      struct Impl {
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void potrf_upper(MemberType &member, 
                         const int m, 
                         T *__restrict__ A, const int as0, const int as1,
                         int *info) {
          if (m <= 0) return;

          typedef ArithTraits<T> arith_traits;
          for (int p=0;p<m;++p) {
            const int jend = m-p-1;
            
            T
              *__restrict__ alpha11 = A+(p  )*as0+(p  )*as1,
              *__restrict__ a12t    = A+(p  )*as0+(p+1)*as1,
              *__restrict__ A22     = A+(p+1)*as0+(p+1)*as1;
            
            Kokkos::single(Kokkos::PerTeam(member), [&] (){
                *alpha11 = sqrt(arith_traits::real(*alpha11));                
              });
            member.team_barrier();
            const auto alpha = arith_traits::real(*alpha11);
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member,jend),[&](const int &j) {
                a12t[j*as1] /= alpha;
              });
            member.team_barrier();
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,jend),[&](const int &j) {
                const T aa = arith_traits::conj(a12t[j*as1]);
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,j+1),[&](const int &i) {
                    const T bb = a12t[i*as1];
                    A22[i*as0+j*as1] -= aa*bb;
                  });
              });
            member.team_barrier();
          }
        }

        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void sytrf_lower(MemberType &member, 
                         const int m, 
                         T *__restrict__ A, const int as0, const int as1,
                         int *__restrict__ fpiv, 
                         T *__restrict__ W,
                         int *info) {
          if (m <= 0) return;
          
          typedef ArithTraits<T> arith_traits;
          for (int p=0;p<m;++p) {
            const int iend = m-p-1;
            
            T
              *__restrict__ alpha11 = A+(p  )*as0+(p  )*as1,
              *__restrict__ a21     = A+(p+1)*as0+(p  )*as1,
              *__restrict__ A22     = A+(p+1)*as0+(p+1)*as1;
            
            const auto alpha = arith_traits::real(*alpha11);
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member,iend),[&](const int &i) {
                a21[i*as0] /= alpha;
              });
            member.team_barrier();
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,iend),[&](const int &i) {
                const T aa = a21[i*as0];
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,i+1),[&](const int &j) {
                    const T bb = a21[j*as0];
                    A22[i*as0+j*as1] -= alpha*aa*bb;
                  });
              });
            member.team_barrier();
          }
        }
      };

      template<typename MemberType>
      static 
      KOKKOS_INLINE_FUNCTION
      void potrf(MemberType &member,
                 const char uplo,
                 const int m, 
                 /* */ T *__restrict__ A, const int lda,
                 int *info) {
        switch (uplo) {
        case 'U':
        case 'u': {
          Impl::potrf_upper(member, 
                            m,
                            A, 1, lda,
                            info); 
          break;
        }
        case 'L':
        case 'l': {
          Kokkos::abort("not implemented");
          break;
        }
        default:
          Kokkos::abort("Invalid uplo character");
        }
      }


      template<typename MemberType>
      static 
      KOKKOS_INLINE_FUNCTION
      void sytrf(MemberType &member,
                 const char uplo,
                 const int m, 
                 /* */ T *__restrict__ A, const int lda,
                 /* */ int *__restrict__ P,
                 /* */ T *__restrict__ W,
                 int *info) {
        switch (uplo) {
        case 'U':
        case 'u': {
          Kokkos::abort("not implemented");
          break;
        }
        case 'L':
        case 'l': {
          Impl::sytrf_lower(member, 
                            m,
                            A, 1, lda,
                            P,
                            W,
                            info); 
          break;
        }
        default:
          Kokkos::abort("Invalid uplo character");
        }
      }

    };

}

#endif
