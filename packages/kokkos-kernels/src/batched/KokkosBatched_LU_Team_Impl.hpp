#ifndef __KOKKOSBATCHED_LU_TEAM_IMPL_HPP__
#define __KOKKOSBATCHED_LU_TEAM_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_LU_Team_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {
    ///
    /// Team Impl
    /// =========

    ///
    /// LU no piv
    ///
    
    template<typename MemberType>
    struct TeamLU<MemberType,Algo::LU::Unblocked> {
      template<typename AViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, const AViewType &A,
             const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0) {
        return TeamLU_Internal<Algo::LU::Unblocked>::invoke(member,
                                                            A.extent(0), A.extent(1),
                                                            A.data(), A.stride_0(), A.stride_1(),
                                                            tiny);
      }
    };
    
    template<typename MemberType>
    struct TeamLU<MemberType,Algo::LU::Blocked> {
      template<typename AViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, const AViewType &A,
             const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0) {
        return TeamLU_Internal<Algo::LU::Blocked>::invoke(member,
                                                          A.extent(0), A.extent(1),
                                                          A.data(), A.stride_0(), A.stride_1(),
                                                          tiny);
      }
    };

  }
}

#endif
