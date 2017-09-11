#ifndef __KOKKOSBATCHED_GEMV_TEAM_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMV_TEAM_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

#include "KokkosBatched_InnerMultipleDotProduct_Serial_Impl.hpp"


namespace KokkosBatched {
  namespace Experimental {

    ///
    /// Team Internal Impl
    /// ====================
    namespace Team {

      template<typename ArgAlgo>
      struct GemvInternal {
        template<typename MemberType,
                 typename ScalarType,
                 typename ValueType>
        KOKKOS_INLINE_FUNCTION
        static int
        invoke(const MemberType &member,
               const int m, const int n, 
               const ScalarType alpha,
               const ValueType *__restrict__ A, const int as0, const int as1,
               const ValueType *__restrict__ x, const int xs0, 
               const ScalarType beta,
               /**/  ValueType *__restrict__ y, const int ys0);
      };

      template<>
      template<typename MemberType,
               typename ScalarType,
               typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int
      GemvInternal<Algo::Gemv::Unblocked>::
      invoke(const MemberType &member,
             const int m, const int n, 
             const ScalarType alpha,
             const ValueType *__restrict__ A, const int as0, const int as1,
             const ValueType *__restrict__ x, const int xs0,
             const ScalarType beta,
             /**/  ValueType *__restrict__ y, const int ys0) {
        // y = beta y + alpha A x
        // y (m), A(m x n), B(n)

        typedef ValueType value_type;

        if      (beta == 0) TeamSetInternal  ::invoke(member, m, 0,    y, ys0);
        else if (beta != 1) TeamScaleInternal::invoke(member, m, beta, y, ys0);
      
        if (alpha != 0) {
          if (m <= 0 || n <= 0) return 0;

          const value_type alpha_value(alpha);
          if (beta != 1) 
            member.barrier();
            
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,m),[&](const int &i) {
              value_type t(0);
              const value_type *__restrict__ tA = (A + i*as0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
              for (int j=0;j<n;++j)
                t += tA[j*as1]*x[j*xs0];
              y[i*ys0] += alpha_value*t;
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
      GemvInternal<Algo::Gemv::Blocked>::
      invoke(const MemberType &member,
             const int m, const int n, 
             const ScalarType alpha,
             const ValueType *__restrict__ A, const int as0, const int as1,
             const ValueType *__restrict__ x, const int xs0,
             const ScalarType beta,
             /**/  ValueType *__restrict__ y, const int ys0) {
        // y = beta y + alpha A x
        // y (m), A(m x n), B(n)

        typedef ValueType value_type;
        enum : int {
          mbAlgo = Algo::Gemv::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>()
        };
          
        if      (beta == 0) TeamSetInternal  ::invoke(member, m, 0,    y, ys0);
        else if (beta != 1) TeamScaleInternal::invoke(member, m, beta, y, ys0);
      
        if (alpha != 0) {
          if (m <= 0 || n <= 0) return 0;

          const value_type alpha_value(alpha);
          if (beta != 1) 
            member.barrier();
        
          InnerMultipleDotProduct<mbAlgo> inner(as0, as1, xs0, ys0);
          const int tsize = member.team_size();
          const int mb = min(m/tsize + (m%tsize)>0, mbAlgo), mp = m%mb;
            
          Kokkos::parallel_for
            (Kokkos::TeamThreadRange(member, (m/mb) + (mp>0)),
             [&](const int &ii) {
              const int i = ii*mb;
              inner.serial_invoke(alpha_value, A+i*as0, x, (i+mb) > m ? (m-i) : mb, n, y+i*ys0 );
            });
          member.team_barrier();
        }
          
        return 0;
      }
    }
  } 
}


#endif
